#DATA READING
branch = CSV.read("branchdata_$system_name.csv"; delim=',')
bus = CSV.read("busdata_$system_name.csv"; delim=',')
# branch=readdlm("branchdata_$system_name.csv",';')
# bus=readdlm("busdata_$system_name.csv",';')
using Printf
using CSV
using LinearAlgebra
#SYSTEM PARAMETERS
# Sbase=bus[1,4]
Sbase = 100
Vbase=1.00
nbus=size(bus,1)-3
nbranch=size(branch,1)-2
# DEFINE PARAMETERS DE NODOS
mutable struct DATA1
 	busnum
 	bustype
 	V0
 	Th0
 	Pd
 	Qd
 	Pg0
 	Qg0
 	Vg
 	gshb
 	bshb
	Qmax
	Qmin
	Pmax
	Pmin
end
# DEFINE PARAMETERS DE LINEAS
mutable struct DATA2
	branchnum
	from
	to
	r
	x
	c
	a
	fi
	slmax
	z2
	g
	b
	tmin
	tmax
end
# #DEFINE PARAMETERS Bus.DATA1  Branch.DATA2
Bus=DATA1(bus[4:nbus+3,1], bus[4:nbus+3,5], bus[4:nbus+3,6], bus[4:nbus+3,7], bus[4:nbus+3,8], bus[4:nbus+3,9],
bus[4:nbus+3,10], bus[4:nbus+3,11],bus[4:nbus+3,13],bus[4:nbus+3,16], bus[4:nbus+3,17], bus[4:nbus+3,14],bus[4:nbus+3,15],bus[4:nbus+3,19],bus[4:nbus+3,20])
Branch=DATA2(collect(1:nbranch),branch[3:nbranch+2,1], branch[3:nbranch+2,2], branch[3:nbranch+2,7],
branch[3:nbranch+2,8], branch[3:nbranch+2,9], branch[3:nbranch+2,15],branch[3:nbranch+2,16],
branch[3:nbranch+2,22],zeros(Float64,1,nbranch), zeros(Float64,1,nbranch),zeros(Float64,1,nbranch),branch[3:nbranch+2,17],branch[3:nbranch+2,18])

#CALCULATE PARAMETERS
for i in 1:3
	# Bus.busLoc[i] = i
	Bus.Pd[i] = Bus.Pd[i]/(Sbase)
	# Bus.Qd[i] /= Sbase
	# Bus.Pg0[i] /= Sbase
	# Bus.Qg0[i] /= Sbase
end
# Bus.Pd ./=(Sbase)
Bus.Qd = Bus.Qd / Sbase
Bus.Pg0 =Bus.Pg0 / Sbase
Bus.Qg0 = Bus.Qg0 / Sbase
Bus.Th0 = Bus.Th0*pi/180
Branch.fi = Branch.fi*pi/180
Branch.slmax = Branch.slmax./Sbase
Branch.c=Branch.c / 2
Branch.z2=Branch.r.^2+Branch.x.^2
Busnumber=(collect(1:nbus))
Bus.Qmax=Bus.Qmax / Sbase
Bus.Qmin=Bus.Qmin / Sbase
Bus.Pmax=Bus.Pmax / Sbase
Bus.Pmin=Bus.Pmin / Sbase

for i=1:nbranch
	Branch.g[i]=Branch.r[i]/(Branch.r[i]^2+Branch.x[i]^2)
	Branch.b[i]=-Branch.x[i]/(Branch.r[i]^2+Branch.x[i]^2)

	if Branch.a[i]==0
	Branch.a[i]=1
	else
	Branch.a[i]=1/Branch.a[i]
	end
end
if IData==1
#MUESTRA DE PARAMETROS
@printf "------------------------------------------------------------------------------\n"
	@printf "---------------------------------INPUT DATA-----------------------------------\n"
	@printf "------------------------------------------------------------------------------\n"
	@printf "                                  BUS_DATA                                    \n"
	@printf "------------------------------------------------------------------------------\n"
	@printf "   Bus  Type    V0       Th0        Pd        Qd        Pg0        Qg0        \n"
	@printf "------------------------------------------------------------------------------\n"
	for i in 1:nbus
		@printf "%5d" float(Busnumber[i])
		@printf "%5d" float(Bus.busnum[i])
		@printf "%5d" float(Bus.bustype[i])
		@printf "%10.4f" float(Bus.V0[i])
		@printf "%10.4f" float(Bus.Th0[i])
		@printf "%10.4f" float(Bus.Pd[i])
		@printf "%10.4f" float(Bus.Qd[i])
		@printf "%10.4f" float(Bus.Pg0[i])
		@printf "%10.4f" float(Bus.Qg0[i])
		@printf "\n"
	end
	@printf "-----------------------------------------------------------------------------------------\n"
	@printf "                                 BRANCH_DATA                                             \n"
	@printf "-----------------------------------------------------------------------------------------\n"
	@printf " Branch From  To     r         x         bshb        a        fi       smax \n"
	@printf "-----------------------------------------------------------------------------------------\n"
	#change of data
	Branchfrom=zeros(Float64,1,nbranch)
	Branchto=zeros(Float64,1,nbranch)
	for i=1:nbranch
	a=find(isequal(Branch.from[i]),Bus.busnum)
	Branchfrom[i]=a[1]
	b=find(isequal(Branch.to[i]),Bus.busnum)
	Branchto[i]=b[1]
	end
	for i in 1:nbranch
		@printf "%5d" float(i)
		@printf "%5d" float(Branch.from[i])
		@printf "%5d" float(Branchfrom[i])
		@printf "%5d" float(Branch.to[i])
		@printf "%5d" float(Branchto[i])
		@printf "%10.4f" float(Branch.r[i])
		@printf "%10.4f" float(Branch.x[i])
		@printf "%10.4f" float(Branch.c[i])
		@printf "%10.4f" float(Branch.a[i])
		@printf "%10.4f" float(Branch.fi[i])
		@printf "%10.4f" float(Branch.slmax[i])
		@printf "\n"
	end
	@printf "-----------------------------------------------------------------------------------------\n"
	@printf "                                                                                         \n"
else
	#change of data
	Branchfrom=zeros(Float64,1,nbranch)
	Branchto=zeros(Float64,1,nbranch)
	for i=1:nbranch
	a=find(isequal(Branch.from[i]),Bus.busnum)
	Branchfrom[i]=a[1]
	b=find(isequal(Branch.to[i]),Bus.busnum)
	Branchto[i]=b[1]
	end
end


#Identificacion de lineas que salen y que llegan a una determinada barra
	in_lines = [Int[] for i in 1:nbus] # indices of incoming branches
	out_lines = [Int[] for i in 1:nbus] # indices of outgoing branches

	for i in 1:nbranch

		a = convert(Int, Branchfrom[i])
		b = convert(Int, Branchto[i])

	    push!(out_lines[a],i)
	    push!(in_lines[b],i)
	end
#location of Transformers
	TrafoLocation = zeros(1,nbranch)
	for i=1:nbranch
		if Branch.a[i]!= 1
	       TrafoLocation[i] =+i
	    else
	    end
	end
#Location of Shunt Compensation
	QShuntLocation = zeros(1,nbranch)
	for i=1:nbus
		if Bus.bshb[i]!= 0
	       QShuntLocation[i] +=1
	    else
	    end
	end
	nQshunt=convert(Int, sum(QShuntLocation))
# Number of Generators
	ngen=0
	GenLocation= zeros(1,nbus)
	for i=1:nbus
	    if Bus.bustype[i] != 0
	        ngen+=1
			GenLocation[i] =+i
	    else
	    end
	end
#Number of Trafos
	nTrafo=0
	for i=1:nbranch
	    if TrafoLocation[i]!=0
	        nTrafo+=1
	    else
	    end
	end
