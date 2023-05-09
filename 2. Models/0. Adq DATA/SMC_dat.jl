#DATA READING
# file_name = "IEEE14"
# Sbase = 100
# IData = 1
using JuMP, MATLAB
#SYSTEM Parameters
# Sbase = 100
# Vnom = 1.00
# Vmax = 1.05
# Vmin = 0.95
# pfcap = 0.95
# pfind = 0.99


# Base System
Sbase = 100
# Show Input Data
IData = 1


system_name = "IEEE14"


using Printf
using CSV
using Graphs, LightGraphs
bus = CSV.read("bus_data_$system_name.csv"; delim=',')
branch = CSV.read("branch_data_$system_name.csv"; delim=',')
cgen = CSV.read("cgen_data_$system_name.csv"; delim=',')
wgen = CSV.read("wgen_data_$system_name.csv"; delim=',')
scen = CSV.read("scen_data_$system_name.csv"; delim=',')
nbus=size(bus,1)
nbranch=size(branch,1)
ncgen=size(cgen,1)
nwgen=size(wgen,1)
nscen=size(scen,1)
@show nbus
@show nbranch
@show ncgen
@show nwgen
@show nscen
# dbc1 = 1 #datos adicionales para la estructura DataBus
dbc2 = 3 #datos adicionales para la estructura DataBranch
# databuscusto = eye(nbus,dbc1) # To add Columns Customized

# databranchcusto = eye(nbranch,dbc2) # To add Columns Customized
# databranchcusto = ones(length(nbranch),length(dbc2)) # To add Columns Customized
databranchcusto = ones(nbranch,dbc2) # To add Columns Customized
@show databranchcusto
#

mutable struct DATA1
	busnum
	bustype
	busname
	V0
	Th0
	Pd
	Qd
	Pg0
	Qg0
	gshb
	bshb
	n0
	Vgen
	scost
	busLoc
end
# struct Foo
#            bar
#            baz::Int
#            qux::Float64
#        end
mutable struct DATA2
	branchnum
	from
	to
	r
	x
	c
	a
	fi
	smax
	branchtype
	reg
	ntmin
	ntmax
	n0
	branchLoc
	g
	b
	z2
end
#
mutable struct DATA3
	cgennum
	bus
	Pg0
	Qg0
	Pg_max
	Pg_min
	Qg_max
	Qg_min
	Vg
	ag
	bg
	cg
	status
	shedding
end
#
mutable struct DATA4
	wgennum
	bus
	Pg0
	Qg0
	Pg_max
	Pg_min
	Qg_max
	Qg_min
	Vg
	ag
	bg
	cg
	status
	scenarios
end
#
mutable struct DATA5
	scennum
	scenID
	wgen
	P_W_RT
	Prob
end
#
Bus=DATA1(bus[:,2], bus[:,6], bus[:,3],	bus[:,7], bus[:,8], bus[:,9], bus[:,10], bus[:,11], bus[:,12], bus[:,17], bus[:,18], bus[:,19], bus[:,14],bus[:,21], bus[:,1])
Branch=DATA2(branch[:,1], branch[:,2], branch[:,3],
	branch[:,8], branch[:,9], branch[:,10], branch[:,16], branch[:,17], branch[:,23],
	branch[:,7], branch[:,12], branch[:,18], branch[:,19], branch[:,6], branch[:,1], databranchcusto[:,1], databranchcusto[:,2], databranchcusto[:,3])

Cgen=DATA3(cgen[:,1], cgen[:,2], cgen[:,3],	cgen[:,4], cgen[:,5], cgen[:,6], cgen[:,7], cgen[:,8], cgen[:,9],
	cgen[:,10], cgen[:,11], cgen[:,12], cgen[:,13], cgen[:,14])
#
Wgen=DATA4(wgen[:,1], wgen[:,2], wgen[:,3],	wgen[:,4], wgen[:,5], wgen[:,6], wgen[:,7], wgen[:,8], wgen[:,9],
		wgen[:,10], wgen[:,11], wgen[:,12], wgen[:,13], wgen[:,14])
Scen=DATA5(scen[:,1], scen[:,2], scen[:,3],	scen[:,4], scen[:,5])
#
# calculated parameters
for i in 1:nbranch
	if Branch.r[i] == 0
		Branch.r[i] = 1e-6
	end
end

for i in 1:nbus
	Bus.busLoc[i] = i
	# Bus.Pd[i] /= Sbase
	# Bus.Qd[i] /= Sbase
	# Bus.Pg0[i] /= Sbase
	# Bus.Qg0[i] /= Sbase
end
Bus.Pd /= Sbase
Bus.Qd /= Sbase
Bus.Pg0 /= Sbase
Bus.Qg0 /= Sbase

# for i in 1:nbranch
# 	Branch.z2[i] = Branch.r[i]^2+Branch.x[i]^2
# 	Branch.g[i] = Branch.r[i]/(Branch.r[i]^2+Branch.x[i]^2)
# 	Branch.b[i] = -Branch.x[i]/(Branch.r[i]^2+Branch.x[i]^2)
# 	Branch.c[i] /= 2
# 	if Branch.a[i] == 0
# 		Branch.a[i] = 1
# 	else
# 		Branch.a[i] = 1/Branch.a[i]
# 	end
# 	Branch.fi[i] *= 3.14159265359/180
# 	Branch.smax[i] = Branch.smax[i]/Sbase
#
# end

for i in 1:nbranch

	Branch.z2[i] = Branch.r[i]^2+Branch.x[i]^2
	Branch.g[i] = Branch.r[i]/(Branch.r[i]^2+Branch.x[i]^2)
	Branch.b[i] = -Branch.x[i]/(Branch.r[i]^2+Branch.x[i]^2)
	Branch.c[i] /= 2
	if Branch.a[i] == 0
		Branch.a[i] = 1
	else
		Branch.a[i] = 1/Branch.a[i]
	end
	Branch.fi[i] *= 3.14159/180
end
Branch.smax /= Sbase

Cgen.Pg0 /= Sbase
Cgen.Qg0 /= Sbase
Cgen.Pg_max /= Sbase
Cgen.Pg_min /= Sbase
Cgen.Qg_max /= Sbase
Cgen.Qg_min /= Sbase

Wgen.Pg0 /= Sbase
Wgen.Qg0 /= Sbase
Wgen.Pg_max /= Sbase
Wgen.Pg_min /= Sbase
Wgen.Qg_max /= Sbase
Wgen.Qg_min /= Sbase

Scen.P_W_RT /= Sbase

for i in 1:nbranch
	Branch.g[i] = Branch.g[i]*Branch.n0[i]
	Branch.b[i] = Branch.b[i]*Branch.n0[i]
	Branch.c[i] = Branch.c[i]*Branch.n0[i]
	Branch.smax[i] = Branch.smax[i]*Branch.n0[i]
end
for i in 1:nbus
	Bus.bshb[i] = Bus.bshb[i]*Bus.n0[i]
end
if IData == 1
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "-------------------------------------------------INPUT-DATA----------------------------------------------\n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "                                                  BUS-DATA                                                     \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "   Bus  Type    V0      Th0      Pd      Qd     Pg0     Qg0     bshb    n0  scost\n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	for i in 1:nbus
		@printf "%5d" float(i)
		@printf "%5d" float(Bus.bustype[i])
		@printf "%10.4f" float(Bus.V0[i])
		@printf "%8.2f" float(Bus.Th0[i])
		@printf "%8.2f" float(Bus.Pd[i])
		@printf "%8.2f" float(Bus.Qd[i])
		@printf "%8.2f" float(Bus.Pg0[i])
		@printf "%8.2f" float(Bus.Qg0[i])
		@printf "%8.2f" float(Bus.bshb[i])
		@printf "%5d" float(Bus.n0[i])
		@printf "%5d" float(Bus.scost[i])
		@printf "\n"
	end
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "                                                 BRANCH-DATA                                                     \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf " Branch From  To     r         x         g        b         bshb      a      smax   n0 \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	#change of data
	Branchfrom=zeros(Float64,1,nbranch)
	Branchto=zeros(Float64,1,nbranch)
	for i=1:nbranch
		a=findall(isequal(Branch.from[i]),Bus.busnum)
		# a=find(isequal(Branch.from[i]),Bus.busnum)
		Branchfrom[i]=a[1]
		b=findall(isequal(Branch.to[i]),Bus.busnum)
		# b=find(isequal(Branch.to[i]),Bus.busnum)
		Branchto[i]=b[1]
	end
	for i in 1:nbranch
		@printf "%5d" float(i)
		@printf "%5d" float(Branch.from[i])
	    @printf "%5d" float(Branch.to[i])
		@printf "%10.4f" float(Branch.r[i])
		@printf "%10.4f" float(Branch.x[i])
		# @printf "%10.4f" float(Branch.z2[i])
		@printf "%10.4f" float(Branch.g[i])
		@printf "%10.4f" float(Branch.b[i])
		@printf "%10.4f" float(Branch.c[i])
		@printf "%8.2f" float(Branch.a[i])
		# @printf "%6.0f" float(Branch.fi[i])
		@printf "%8.2f" float(Branch.smax[i])
		@printf "%4.0f" float(Branch.n0[i])
		@printf "\n"
	end
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "                                                                                                         \n"
	@printf "                                                  CGEN-DATA                                                     \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "    N°  Bus   Pg0     Qg0       Pg_max    Pg_min    Qg_max    Qg_min   ag   bg    cg  Stat  shed\n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	for i in 1:ncgen
		@printf "%5d" float(i)
	    @printf "%5d" float(Cgen.bus[i])
		@printf "%8.2f" float(Cgen.Pg0[i])
		@printf "%10.2f" float(Cgen.Qg0[i])
		@printf "%10.2f" float(Cgen.Pg_max[i])
		@printf "%10.2f" float(Cgen.Pg_min[i])
		@printf "%10.2f" float(Cgen.Qg_max[i])
		@printf "%10.2f" float(Cgen.Qg_min[i])
		@printf "%5d" float(Cgen.ag[i])
		@printf "%5d" float(Cgen.bg[i])
		@printf "%5d" float(Cgen.cg[i])
		@printf "%5d" float(Cgen.status[i])
		@printf "%5d" float(Cgen.shedding[i])
		@printf "\n"
	end
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "                                                                                                         \n"
	@printf "                                                  WGEN-DATA                                                     \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "    N°  Bus   Pg0     Qg0       Pg_max    Pg_min    Qg_max    Qg_min      bg   scena  status\n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	for i in 1:nwgen
		@printf "%5d" float(i)
	    @printf "%5d" float(Wgen.bus[i])
		@printf "%8.4f" float(Wgen.Pg0[i])
		@printf "%10.4f" float(Wgen.Qg0[i])
		@printf "%10.4f" float(Wgen.Pg_max[i])
		@printf "%10.4f" float(Wgen.Pg_min[i])
		@printf "%10.4f" float(Wgen.Qg_max[i])
		@printf "%10.4f" float(Wgen.Qg_min[i])
		@printf "%8d" float(Wgen.bg[i])
		@printf "%5d" float(Wgen.scenarios[i])
		@printf "%8d" float(Wgen.status[i])
		@printf "\n"
	end
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "                                                  SCEN-DATA                                                     \n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	@printf "    N°  Scen Wgen   P_W      Prob\n"
	@printf "---------------------------------------------------------------------------------------------------------\n"
	for i in 1:nscen
		@printf "%5d" float(i)
	    @printf "%5d" float(Scen.scenID[i])
		@printf "%5d" float(Scen.wgen[i])
		@printf "%10.4f" float(Scen.P_W_RT[i])
		@printf "%10.4f" float(Scen.Prob[i])
		@printf "\n"
	end
	@printf "---------------------------------------------------------------------------------------------------------\n"
else
	#change of data
	Branchfrom=zeros(Float64,1,nbranch)
	Branchto=zeros(Float64,1,nbranch)
	for i=1:nbranch
	a=findall(isequal(Branch.from[i]),Bus.busnum)
	Branchfrom[i]=a[1]
	b=findall(isequal(Branch.to[i]),Bus.busnum)
	Branchto[i]=b[1]
	end
end

in_lines = [Int[] for i in 1:nbus] # indices of incoming branches
out_lines = [Int[] for i in 1:nbus] # indices of outgoing branches

for i in 1:nbranch

	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])

    push!(out_lines[a],i)
    push!(in_lines[b],i)

end

@show in_lines
@show out_lines

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
@show Y
# @mput Y

# mat"A=real($Y) + speye(size($Y))"
mat"A=abs(imag($Y)) + speye(size($Y))"
mat"p=amd(A)"
mat"L = chol(A,'lower')"
mat"Ly = chol(real($Y),'lower')"
mat"Lp = chol(A(p,p),'lower')"
mat"Ap = A(p,p)"
mat"figure"
mat"hold on"  # evaluate a MATLAB function
mat"subplot(2,3,1)"
mat"spy($Y)"
mat"title('Original Matrix Y')"
mat"subplot(2,3,2)"
mat"spy(A)"
mat"title('Matrix A')"
mat"subplot(2,3,3)"
mat"spy(A(p,p))"
mat"title('AMD ordered Y')"
mat"subplot(2,3,4)"
mat"spy(Ly)"
mat"title('Cholesky factor of Y')"
mat"subplot(2,3,5)"
mat"spy(L)"
mat"title('Cholesky factor of A')"
mat"subplot(2,3,6)"
mat"spy(Lp)"
mat"title('Cholesky factor of AMD ordered A')"
mat"grid on"  # evaluate a MATLAB function
mat"hold off"  # evaluate a MATLAB function
mat"figure"
mat"hold on"  # evaluate a MATLAB function
mat"spy(Lp*Lp')"
mat"hold off"  # evaluate a MATLAB function
# y = range(2.0, stop=3.0, length=500)
# mat"""
#     $u = $x + $y
# 	$v = $x - $y
# """
@mget Ap
@mget A
@mget Lp
@show A               # u and v are accessible from Julia
@show Ap               # u and v are accessible from Julia
# As = Dict()
row, col = size(Lp)
As = zeros(row,col)
for i in 1:col
	for j in 1:row
		if Lp[i,j] != 0
			As[i,j] = 1
		end
	end
end
@show As
mat"figure"
mat"hold on"  # evaluate a MATLAB function
mat"spy($As)"
mat"hold off"  # evaluate a MATLAB function
# # g = Graph(10,20)
# Am = adjacency_matrix(Ap)
# @show Am
# matrixM = maximal_cliques(Am)
A  = As*As'
include("MaximalClique.jl")
