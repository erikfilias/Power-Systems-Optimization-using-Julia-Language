# =============================================================================#
#Description: File to read Data
#Date: 08.14.2019
#Version: 0.0.1.081419
#Model: Non Lineal Rectangular
# =============================================================================#
println("***  Start Data Reading ***")
bus    = CSV.read("bus_data_$file_name.csv")
print(bus)
# bus    = CSV.read("bus_data_$file_name.csv"; delim=',')
branch = CSV.read("branch_data_$file_name.csv")
# branch = CSV.read("branch_data_$file_name.csv"; delim=',')
nbus=size(bus,1)
nbranch=size(branch,1)
#dbc1 = 1 #datos adicionales para la estructura DataBus
dbc2 = 3 #datos adicionales para la estructura DataBranch
#databuscusto = zeros(nbus,dbc1) # To add Columns Customized
databranchcusto = zeros(nbranch,dbc2) # To add Columns Customized


# Structure for Bus
mutable struct DATA1
	busId
	bustype
	Vgen
	Th0
	Pd
	Qd
	Pg0
	Qg0
	gshb
	bshb
	busLoc # Location or Ubiety of the branch on the matrix
end

# Structure for Line
mutable struct DATA2
		from
		to
		r
		x
		c
		a
		fi
		branchLoc  # Location or Ubiety of the branch on the matrix
		g
		b
		z2
end

Bus=DATA1(
	        bus[:,2],
			bus[:,3],
			bus[:,4],
			bus[:,5],
			bus[:,6],
			bus[:,7],
			bus[:,8],
			bus[:,9],
			bus[:,10],
			bus[:,11],
			bus[:,1]) # Location or Ubiety of the bus on the matrix


Branch=DATA2(
		 branch[:,2],
		 branch[:,3],
		 branch[:,4],
		 branch[:,5],
		 branch[:,6],
		 branch[:,7],
		 branch[:,8],
		 branch[:,1], # Location or Ubiety of the branch on the matrix
		 databranchcusto[:,1], # calculation of G from line
		 databranchcusto[:,2], # calculation of B from line
		 databranchcusto[:,3])  # calculation of Z2 from line


# Parameters for denoted in pu
Sbase = 100
Vmin = 0.95; # magnitude de tensão mínima em pu
Vmax = 1.05; # magnitude de tensão máxima em pu
Vnom = 1.00; # magnitude de tensão nominal em pu

# Modifying the data of Bus
# Power Active & Reactive denoted in pu

Bus.Th0 *=  3.14159265359/180
Bus.Pd /= Sbase
Bus.Qd /= Sbase
Bus.Pg0 /= Sbase
Bus.Qg0 /= Sbase
Bus.busLoc = collect(1:nbus)


# Modifying the data of Line
# Considering the values of fi in degrees
Branch.branchLoc = collect(1:nbranch)
#Branch.c  /=  2
#Branch.a  +=  zeros(nbranch,1)
# Branch.fi *=  3.14159265359/180
# Branch.z2 +=  (Branch.r.^2 + Branch.x.^2)
# Branch.g   =  (Branch.r)./(Branch.z2)
# Branch.b   = -(Branch.x)./(Branch.z2)

for i in 1:nbranch
	Branch.z2[i] = Branch.r[i]^2+Branch.x[i]^2
	Branch.g[i] = Branch.r[i]/(Branch.r[i]^2+Branch.x[i]^2)
	Branch.b[i] = -Branch.x[i]/(Branch.r[i]^2+Branch.x[i]^2)
	Branch.c[i] /= 2
	if Branch.a[i] == 0; Branch.a[i] = 1
	else Branch.a[i] = 1/Branch.a[i]
	end
	Branch.fi[i] *= 3.14159/180
end

# ubicando las lines de salida y entradas de las barras


Branchfrom=zeros(Float64,1,nbranch)
Branchto=zeros(Float64,1,nbranch)
for i=1:nbranch
	a=findall(isequal(Branch.from[i]),Bus.busId)
	Branchfrom[i]=a[1]
	b=findall(isequal(Branch.to[i]),Bus.busId)
	Branchto[i]=b[1]
end


in_lines = [Int[] for i in 1:nbus] # indices of incoming branches
out_lines = [Int[] for i in 1:nbus] # indices of outgoing branches

for i in 1:nbranch

	a = convert(Int, Branchfrom[i])
	b = convert(Int, Branchto[i])

    push!(out_lines[a],i)
    push!(in_lines[b],i)

end
println("***  Finsh Data Reading ***")
