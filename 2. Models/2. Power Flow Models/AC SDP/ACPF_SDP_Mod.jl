# =============================================================================#
#Description: Model SDP
#Date: 08.14.2019
#Version: 0.0.1.081419
#Model: AC SDP
# =============================================================================#


#  ============= Declaration of Variables to Bus ==============================#
@variable(model, x5[1:2*nbus,1:2*nbus],PSD)	# Group of buses:
@variable(model, Pg[Bus.busLoc])
@variable(model, Qg[Bus.busLoc])
@variable(model, Qshunt[Bus.busLoc])
for i in 1:nbus
   set_start_value(Pg[i],Bus.Pg0[i])	# Power Active Supply on node i
   set_start_value(Qg[i],Bus.Qg0[i])	# Power Reactive Supply on node i
   set_start_value(Qshunt[i], 0)		# Aux Var to Reactive Compensation
end
#  ===========  Declaration of Variables to Branch ============================#
@variable(model, Pde[Branch.branchLoc])
@variable(model, Ppa[Branch.branchLoc])
@variable(model, Qde[Branch.branchLoc])
@variable(model, Qpa[Branch.branchLoc])
#  =========== Initial Values to Variables Branch =============================#
for i in 1:nbranch
	set_start_value(Pde[i], 0)	# Power flow ative from i to j nodes
	set_start_value(Ppa[i], 0)	# Power flow active from j to i nodes
	set_start_value(Qde[i], 0)	# Power flow reative from i to j nodes
	set_start_value(Qpa[i], 0)	# Power flow reactive from j to i nodes
end

println("Declared SD Variables : Ok\n")

#  ============= SD constraint==================================
# @SDconstraint(model, x5>= 0)

#  ============= Objective Function=============================
if FO == 0
	@objective(model,Min, sum(Pg[i] for i=1:nbus if Bus.bustype[i]==3))
elseif FO == 1
	@objective(model,Min, sum(Pg[i] for i=1:nbus if Bus.bustype[i]!=0))
elseif FO == 2
	@objective(model,Min, sum(Pde[i]+Ppa[i] for i=1:nbranch))
end
#  ============= Power Balance ================================
for k=1:nbus
    if Bus.bshb[k]!=0
		auxshunt=zeros(2*nbus,2*nbus)
	    auxshunt[k+k-1,k+k-1]   =Bus.bshb[k]
	    auxshunt[k+k,k+k]       =Bus.bshb[k]
		@constraint(model,Qshunt[k]==dot(auxshunt,x5))

    else
        @constraint(model,Qshunt[k]==0)
    end
end

# Creation Temporal matrix
mtemporal=  zeros(2*nbus,2*nbus)
mtemporal1= zeros(2*nbus,2*nbus)

mtemporal2= zeros(2*nbus,2*nbus)
mtemporal3= zeros(2*nbus,2*nbus)

WPpa =zeros(2*nbus,2*nbus)
WQpa =zeros(2*nbus,2*nbus)

WPde   =zeros(2*nbus,2*nbus)
WQde   =zeros(2*nbus,2*nbus)

for i=1:nbranch
    mtemporal=  zeros(2*nbus,2*nbus)
    mtemporal1= zeros(2*nbus,2*nbus)

    mtemporal2= zeros(2*nbus,2*nbus)
    mtemporal3= zeros(2*nbus,2*nbus)

    WPde   =zeros(2*nbus,2*nbus)
    WQde   =zeros(2*nbus,2*nbus)

    a = convert(Int, Branchfrom[i])
    b = convert(Int, Branchto[i])

	mtemporal[a+a-1,a+a-1]= 2*-(Branch.g[i] * (Branch.a[i]^2))
	mtemporal[a+a,a+a]=     2*-(Branch.g[i] * (Branch.a[i]^2))

	mtemporal2[a+a-1,a+a-1]= 2*(-(Branch.b[i] + Branch.c[i]) * (Branch.a[i]^2))
	mtemporal2[a+a,a+a]=     2*(-(Branch.b[i] + Branch.c[i]) * (Branch.a[i]^2))

        mtemporal1[a+a-1,b+b-1]= Branch.a[i] * Branch.g[i]
        mtemporal1[a+a,b+b]=     Branch.a[i] * Branch.g[i]
        mtemporal1[a+a-1,b+b]=  -Branch.a[i] * Branch.b[i]
        mtemporal1[a+a,b+b-1]=   Branch.a[i] * Branch.b[i]

        mtemporal1[b+b-1,a+a-1]= Branch.a[i] * Branch.g[i]
        mtemporal1[b+b,a+a]=     Branch.a[i] * Branch.g[i]
        mtemporal1[b+b,a+a-1]=  -Branch.a[i] * Branch.b[i]
        mtemporal1[b+b-1,a+a]=   Branch.a[i] * Branch.b[i]

        mtemporal3[a+a-1,b+b-1]= Branch.a[i] * Branch.b[i]
        mtemporal3[a+a,b+b]=     Branch.a[i] * Branch.b[i]
        mtemporal3[a+a-1,b+b]=   Branch.a[i] * Branch.g[i]
        mtemporal3[a+a,b+b-1]=  -Branch.a[i] * Branch.g[i]

        mtemporal3[b+b-1,a+a-1]= Branch.a[i] * Branch.b[i]
        mtemporal3[b+b,a+a]=     Branch.a[i] * Branch.b[i]
        mtemporal3[b+b,a+a-1]=   Branch.a[i] * Branch.g[i]
        mtemporal3[b+b-1,a+a]=  -Branch.a[i] * Branch.g[i]

        WPde=mtemporal+mtemporal1
        WQde=mtemporal2+mtemporal3

	@constraint(model,dot((-0.5)*WPde,x5)==Pde[i])
	@constraint(model,dot((0.5)*WQde,x5)==Qde[i])
end

for i=1:nbranch
    mtemporal=  zeros(2*nbus,2*nbus)
    mtemporal1= zeros(2*nbus,2*nbus)

    mtemporal2= zeros(2*nbus,2*nbus)
    mtemporal3= zeros(2*nbus,2*nbus)

    WPpa   =zeros(2*nbus,2*nbus)
    WQpa   =zeros(2*nbus,2*nbus)

    a = convert(Int, Branchfrom[i])
    b = convert(Int, Branchto[i])

	mtemporal[b+b-1,b+b-1]= 2*-Branch.g[i]
	mtemporal[b+b,b+b]=     2*-Branch.g[i]

	mtemporal2[b+b-1,b+b-1]= 2*-(Branch.b[i] + Branch.c[i])
	mtemporal2[b+b,b+b]=     2*-(Branch.b[i] + Branch.c[i])

        mtemporal1[a+a-1,b+b-1]=Branch.a[i]*Branch.g[i]
        mtemporal1[a+a,b+b]=    Branch.a[i]*Branch.g[i]
        mtemporal1[a+a-1,b+b]=  Branch.a[i]*Branch.b[i]
        mtemporal1[a+a,b+b-1]= -Branch.a[i]*Branch.b[i]

        mtemporal1[b+b-1,a+a-1]=Branch.a[i]*Branch.g[i]
        mtemporal1[b+b,a+a]=    Branch.a[i]*Branch.g[i]
        mtemporal1[b+b,a+a-1]=  Branch.a[i]*Branch.b[i]
        mtemporal1[b+b-1,a+a]= -Branch.a[i]*Branch.b[i]

        mtemporal3[a+a-1,b+b-1]=Branch.a[i]*Branch.b[i]
        mtemporal3[a+a,b+b]=    Branch.a[i]*Branch.b[i]
        mtemporal3[a+a-1,b+b]= -Branch.a[i]*Branch.g[i]
        mtemporal3[a+a,b+b-1]=  Branch.a[i]*Branch.g[i]

        mtemporal3[b+b-1,a+a-1]=Branch.a[i]*Branch.b[i]
        mtemporal3[b+b,a+a]=    Branch.a[i]*Branch.b[i]
        mtemporal3[b+b,a+a-1]= -Branch.a[i]*Branch.g[i]
        mtemporal3[b+b-1,a+a]=  Branch.a[i]*Branch.g[i]

        WPpa=mtemporal+mtemporal1
        WQpa=mtemporal2+mtemporal3

	@constraint(model,dot((-0.5)*WPpa,x5)==Ppa[i])
	@constraint(model,dot((0.5)*WQpa,x5)==Qpa[i])
end

println("Declared Branch Flow constraints : Ok\n")


for k in 1:nbus

    if Bus.bustype[k] !=    0
 		@constraint(model, Pg[k]==sum(Pde[i] for i in out_lines[k])
                   + sum(Ppa[i] for i in in_lines[k])  + Bus.Pd[k])
 		@constraint(model, Qg[k]==-Qshunt[k]
                   + sum(Qde[i] for i in out_lines[k])
                   + sum(Qpa[i] for i in in_lines[k])  +  Bus.Qd[k])
          else

		@constraint(model,sum(Pde[i] for i in out_lines[k])
                  + sum(Ppa[i] for i in in_lines[k]) + Bus.Pd[k] == 0)
		@constraint(model,- Qshunt[k]
                  + sum(Qde[i] for i in out_lines[k])
                  + sum(Qpa[i] for i in in_lines[k]) + Bus.Qd[k] == 0)
    end
end
#  ============= Variables Fixed ================================
for i=1:nbus
    if Bus.bustype[i] ==2
        @constraint(model,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vgen[i]^2)
		@constraint(model,Pg[i]==Bus.Pg0[i])
    end

    if Bus.bustype[i] == 0
        @constraint(model, Qg[i] == 0)
		@constraint(model, Pg[i] == 0)
    end

    if Bus.bustype[i] ==3
		@constraint(model,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vgen[i]^2)
        @constraint(model,x5[i+i,i+i]==0)
    end

end
println("Declared fixed variables : Ok \n")
