#model
m=Model(solver=solver)
#-----------------------------------------------------------------------------
#SD Variables
# Group of active resources
@variable(m, x1[1:2*ngen,1:2*ngen],SDP)
# Group of slack variables for active resources:
@variable(m, x2[1:2*ngen,1:2*ngen],SDP)
# Group of reactive resources:
@variable(m, x3[1:2*ngen,1:2*ngen],SDP)
# Group of slack variables for reactive resources:
@variable(m, x4[1:2*ngen,1:2*ngen],SDP)
# Group of buses:
@variable(m, x5[1:2*nbus,1:2*nbus],SDP)
# Group of slack variables for buses:
@variable(m, x6[1:2*nbus,1:2*nbus],SDP)

println("Declared SD variables : Ok \n")
#------------------------------------------------------------------------------
@SDconstraint(m, x1>= 0)
@SDconstraint(m, x2>= 0)
@SDconstraint(m, x3>= 0)
@SDconstraint(m, x4>= 0)
@SDconstraint(m, x5>= 0)
@SDconstraint(m, x6>= 0)

println("Declared SD constraint : Ok \n")
#-------------------------------------------------------------------------------
#Extra Data

# cte=[0.04302926	20	0;
# 0.25	40	0;
# 0.01	40	0;
# 0.01	40	0;
# 0.01	40	0]
#
# #Datos para 14 Barras
# Vmin=[0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95]
# Vmax=[1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10 1.10]
# Pmax=[2.5 0.4 0 0 0 0 0 0 0 0 0 0 0 0]
# Pmin=[0 0.4 0 0 0 0 0 0 0 0 0 0 0 0]
# Qmax=[0 1 1 0 0 1 0 1 0 0 0 0 0 0]
# Qmin=[-1 0 0 0 0 0 0 0 0 0 0 0 0 0]


# Wg0=zeros(2*ngen,2*ngen)
# for i=1:ngen
#     Wg0[i+i-1,i+i-1]=cte[i,1]
#     Wg0[i+i,i+i]=cte[i,3]
#     Wg0[i+i-1,i+i]=0.5*cte[i,2]
#     Wg0[i+i,i+i-1]=0.5*cte[i,2]
# end

# @objective(m,Min, vecdot(Wg0,x1))#Objective function for OPF
# ------------------------------------------------------------------------------
@objective(m,Min,x1[1,2]) #Objective Function of Power Flow
println("Declared Objetive Function : Ok \n")

# Creation of Ybus to create a funtion Ybus
Branchz=Branch.r+Branch.x*im
Branchy=1 ./Branchz
Y=zeros(nbus,nbus)
Y=complex(Y)
# Elements outside the diagonal of Ybus
for k=1:nbranch
    Y[convert(Int, Branchfrom[k]),convert(Int, Branchto[k])] = Y[convert(Int, Branchfrom[k]),convert(Int, Branchto[k])]+Branchy[k]*Branch.a[k]
    Y[convert(Int, Branchto[k]),convert(Int, Branchfrom[k])] = Y[convert(Int, Branchfrom[k]),convert(Int, Branchto[k])]
end
# Elements inside the diagonal of Ybus
for i=1:nbus
    for j=1:nbranch
        if Branchfrom[j]==i
            Y[i,i] = Y[i,i] + Branchy[j] * Branch.a[j]^2 + Branch.c[j]*im
        elseif Branchto[j]==i
            Y[i,i] = Y[i,i] + Branchy[j] + Branch.c[j]*im
        end
     end
end
for i=1:nbus
        Y[i,i] = Y[i,i] + Bus.bshb[i]*im
end
println("Declared Y Bus : Ok \n")
# ---------------------------------------------------------------------------------------------------------------------------------------
# 1.Power flow equations:
# Creation Temporal matrix
mtemporal=  zeros(2*nbus,2*nbus)
mtemporal1= zeros(2*nbus,2*nbus)

mtemporal2= zeros(2*nbus,2*nbus)
mtemporal3= zeros(2*nbus,2*nbus)

WGp= zeros(2*ngen,2*ngen)
WGq= zeros(2*ngen,2*ngen)

WBp= zeros(2*nbus,2*nbus)
WBq= zeros(2*nbus,2*nbus)

generador=0

for k=1:nbus
    if Bus.bustype[k] !=0

        mtemporal=  zeros(2*nbus,2*nbus)
        mtemporal1= zeros(2*nbus,2*nbus)

        mtemporal2= zeros(2*nbus,2*nbus)
        mtemporal3= zeros(2*nbus,2*nbus)

        WGp= zeros(2*ngen,2*ngen)
        WGq= zeros(2*ngen,2*ngen)

        WBp= zeros(2*nbus,2*nbus)
        WBq= zeros(2*nbus,2*nbus)

        mtemporal[k+k-1,k+k-1]= -2*real(Y[k,k])
        mtemporal[k+k,k+k]=     -2*real(Y[k,k])

        mtemporal2[k+k-1,k+k-1]= 2*imag(Y[k,k])
        mtemporal2[k+k,k+k]=     2*imag(Y[k,k])

        generador+=1

        WGp[generador+generador-1,generador+generador]=1
        WGp[generador+generador,generador+generador-1]=1

        WGq[generador+generador-1,generador+generador]=1
        WGq[generador+generador,generador+generador-1]=1

        a=convert(Int,(size(out_lines[k]))[1])

        for i=1:a
            j=convert(Int,(Branchto[out_lines[k]])[i])

            mtemporal1[k+k-1,j+j-1]=real(Y[k,j])
            mtemporal1[k+k,j+j]=    real(Y[k,j])
            mtemporal1[j+j-1,k+k-1]=real(Y[k,j])
            mtemporal1[j+j,k+k]=    real(Y[k,j])

            mtemporal1[k+k-1,j+j]= -imag(Y[k,j])
            mtemporal1[k+k,j+j-1]=  imag(Y[k,j])
            mtemporal1[j+j,k+k-1]= -imag(Y[k,j])
            mtemporal1[j+j-1,k+k]=  imag(Y[k,j])

            mtemporal3[k+k-1,j+j-1]=-imag(Y[k,j])
            mtemporal3[k+k,j+j]=    -imag(Y[k,j])
            mtemporal3[j+j-1,k+k-1]=-imag(Y[k,j])
            mtemporal3[j+j,k+k]=    -imag(Y[k,j])

            mtemporal3[k+k-1,j+j]=  -real(Y[k,j])
            mtemporal3[k+k,j+j-1]=  real(Y[k,j])
            mtemporal3[j+j,k+k-1]=  -real(Y[k,j])
            mtemporal3[j+j-1,k+k]=  real(Y[k,j])

        end

        b=convert(Int,(size(in_lines[k]))[1])

        for i=1:b
            j=convert(Int,(Branchfrom[in_lines[k]])[i])

            mtemporal1[k+k-1,j+j-1]=real(Y[k,j])
            mtemporal1[k+k,j+j]=    real(Y[k,j])
            mtemporal1[j+j-1,k+k-1]=real(Y[k,j])
            mtemporal1[j+j,k+k]=    real(Y[k,j])

            mtemporal1[k+k-1,j+j]= -imag(Y[k,j])
            mtemporal1[k+k,j+j-1]=  imag(Y[k,j])
            mtemporal1[j+j,k+k-1]= -imag(Y[k,j])
            mtemporal1[j+j-1,k+k]=  imag(Y[k,j])

            mtemporal3[k+k-1,j+j-1]=-imag(Y[k,j])
            mtemporal3[k+k,j+j]=    -imag(Y[k,j])
            mtemporal3[j+j-1,k+k-1]=-imag(Y[k,j])
            mtemporal3[j+j,k+k]=    -imag(Y[k,j])

            mtemporal3[k+k-1,j+j]=  -real(Y[k,j])
            mtemporal3[k+k,j+j-1]=  real(Y[k,j])
            mtemporal3[j+j,k+k-1]=  -real(Y[k,j])
            mtemporal3[j+j-1,k+k]=  real(Y[k,j])

        end
        WBp=mtemporal1+mtemporal
        WBq=mtemporal2+mtemporal3
            @constraint(m,vecdot(0.5*WGp,x1)+vecdot((0.5)*WBp,x5)==Bus.Pd[k])
            @constraint(m,vecdot(0.5*WGq,x3)+vecdot((0.5)*WBq,x5)==Bus.Qd[k])

    elseif Bus.bustype[k] ==0
        mtemporal=  zeros(2*nbus,2*nbus)
        mtemporal1= zeros(2*nbus,2*nbus)

        mtemporal2= zeros(2*nbus,2*nbus)
        mtemporal3= zeros(2*nbus,2*nbus)

        WBp= zeros(2*nbus,2*nbus)
        WBq= zeros(2*nbus,2*nbus)

        mtemporal[k+k-1,k+k-1]= -2*real(Y[k,k])
        mtemporal[k+k,k+k]=     -2*real(Y[k,k])

        mtemporal2[k+k-1,k+k-1]= 2*imag(Y[k,k])
        mtemporal2[k+k,k+k]=     2*imag(Y[k,k])

        a=convert(Int,(size(out_lines[k]))[1])

        for i=1:a
            j=convert(Int,(Branchto[out_lines[k]])[i])

            mtemporal1[k+k-1,j+j-1]=real(Y[k,j])
            mtemporal1[k+k,j+j]=    real(Y[k,j])
            mtemporal1[j+j-1,k+k-1]=real(Y[k,j])
            mtemporal1[j+j,k+k]=    real(Y[k,j])

            mtemporal1[k+k-1,j+j]= -imag(Y[k,j])
            mtemporal1[k+k,j+j-1]=  imag(Y[k,j])
            mtemporal1[j+j,k+k-1]= -imag(Y[k,j])
            mtemporal1[j+j-1,k+k]=  imag(Y[k,j])

            mtemporal3[k+k-1,j+j-1]=-imag(Y[k,j])
            mtemporal3[k+k,j+j]=    -imag(Y[k,j])
            mtemporal3[j+j-1,k+k-1]=-imag(Y[k,j])
            mtemporal3[j+j,k+k]=    -imag(Y[k,j])

            mtemporal3[k+k-1,j+j]=  -real(Y[k,j])
            mtemporal3[k+k,j+j-1]=  real(Y[k,j])
            mtemporal3[j+j,k+k-1]=  -real(Y[k,j])
            mtemporal3[j+j-1,k+k]=  real(Y[k,j])

        end

        b=convert(Int,(size(in_lines[k]))[1])

        for i=1:b
            j=convert(Int,(Branchfrom[in_lines[k]])[i])

            mtemporal1[k+k-1,j+j-1]=real(Y[k,j])
            mtemporal1[k+k,j+j]=    real(Y[k,j])
            mtemporal1[j+j-1,k+k-1]=real(Y[k,j])
            mtemporal1[j+j,k+k]=    real(Y[k,j])

            mtemporal1[k+k-1,j+j]= -imag(Y[k,j])
            mtemporal1[k+k,j+j-1]=  imag(Y[k,j])
            mtemporal1[j+j,k+k-1]= -imag(Y[k,j])
            mtemporal1[j+j-1,k+k]=  imag(Y[k,j])

            mtemporal3[k+k-1,j+j-1]=-imag(Y[k,j])
            mtemporal3[k+k,j+j]=    -imag(Y[k,j])
            mtemporal3[j+j-1,k+k-1]=-imag(Y[k,j])
            mtemporal3[j+j,k+k]=    -imag(Y[k,j])

            mtemporal3[k+k-1,j+j]=  -real(Y[k,j])
            mtemporal3[k+k,j+j-1]=  real(Y[k,j])
            mtemporal3[j+j,k+k-1]=  -real(Y[k,j])
            mtemporal3[j+j-1,k+k]=  real(Y[k,j])

        end
        WBp=mtemporal1+mtemporal
        WBq=mtemporal2+mtemporal3
            @constraint(m,vecdot((0.5)*WBp,x5)==Bus.Pd[k])
            @constraint(m,vecdot((0.5)*WBq,x5)==Bus.Qd[k])
    end
end
# 2.Constraint of reference bus:
for i=1:nbus
    if Bus.bustype[i] !=0
        @constraint(m,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vg[i]^2)
    end
    if Bus.bustype[i] ==3
        @constraint(m,x5[i+i,i+i]==0)

    end
    if Bus.bustype[i] ==2
        @constraint(m,x5[i+i-1,i+i-1]+x5[i+i,i+i]==Bus.Vg[i]*Bus.Vg[i])
    end

end
# 3.Limits of active and reactive power:
g=0
for k=1:nbus
    if Bus.bustype[k] !=0
        g+=1
    @constraint(m,x1[g+g-1,g+g]+x2[g+g-1,g+g-1]==Bus.Pmax[k])
    @constraint(m,x1[g+g-1,g+g]-x2[g+g,g+g]==Bus.Pmin[k])
    @constraint(m,x3[g+g-1,g+g]+x4[g+g-1,g+g-1]==Bus.Qmax[k])
    @constraint(m,x3[g+g-1,g+g]-x4[g+g,g+g]==Bus.Qmin[k])
    # @constraint(m,x1[g+g-1,g+g]+x2[g+g-1,g+g-1]==Pmax[k])
    # @constraint(m,x1[g+g-1,g+g]-x2[g+g,g+g]==Pmin[k])
    # @constraint(m,x3[g+g-1,g+g]+x4[g+g-1,g+g-1]==Qmax[k])
    # @constraint(m,x3[g+g-1,g+g]-x4[g+g,g+g]==Qmin[k])
    end
    if Bus.bustype[k] ==2
    @constraint(m,x1[g+g-1,g+g-1]==Bus.Pg0[k]^2)
    end
end
# 4.Limits of amplitude at each bus:
WV=zeros(2*nbus,2*nbus)
Wf=zeros(2*nbus,2*nbus)
Winf=zeros(2*nbus,2*nbus)
Wsup=zeros(2*nbus,2*nbus)

for k=1:nbus

    WV=zeros(2*nbus,2*nbus)
    Winf=zeros(2*nbus,2*nbus)
    Wsup=zeros(2*nbus,2*nbus)

    WV[k+k-1,k+k-1]=1
    WV[k+k,k+k]=1
    Winf[k+k,k+k]=1
    Wsup[k+k-1,k+k-1]=1

    @constraint(m,vecdot(WV,x5)+vecdot(Wsup,x6)==1.1^2)
    @constraint(m,vecdot(WV,x5)-vecdot(Winf,x6)==0.9^2)
end
# 5. Auxiliary variables
Wauxg=zeros(2*ngen,2*ngen)
Wauxr=zeros(2*ngen,2*ngen)
for k=1:ngen
    Wauxg=zeros(2*ngen,2*ngen)
    Wauxr=zeros(2*ngen,2*ngen)
        Wauxg[k+k,k+k]=1
        Wauxr[k+k,k+k]=1

        @constraint(m,vecdot(Wauxg,x1)==1)
        @constraint(m,vecdot(Wauxr,x3)==1)
end
print(m)
status=solve(m)
