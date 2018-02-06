#DATA READING
#bus = readdlm("IEEE30a.bus")
bus = readdlm("$system_name.bus")
#branch = readdlm("IEEE30a.branch")
branch = readdlm("$system_name.branch")

#SYSTEM Parameters
Sbase = 100
Vnom = 1.00


nbus=size(bus,1)
nbranch=size(branch,1)

type DATA1
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
end

type DATA2
	branchnum
	from
	to
	r
	x
	c
	a
	fi
	smax
	z2
end

#Parameters

Bus=DATA1(bus[:,1], bus[:,2], bus[:,3], bus[:,4], bus[:,5], bus[:,6], bus[:,7], bus[:,8], bus[:,9], bus[:,10], bus[:,11] )
Branch=DATA2(branch[:,1], branch[:,2], branch[:,3], branch[:,4], branch[:,5], branch[:,6], branch[:,7], branch[:,8], branch[:,9], branch[:,1] )

#calculated parameters

for i in 1:nbus
	Bus.Pd[i] /= Sbase
	Bus.Qd[i] /= Sbase
	Bus.Pg0[i] /= Sbase
	Bus.Qg0[i] /= Sbase
end

#branch_z2 = branch_num
for i in 1:nbranch
	Branch.z2[i] = Branch.r[i]^2+Branch.x[i]^2
	Branch.c[i] /= 2
	if Branch.a[i] == 0
		Branch.a[i] = 1
	else
		Branch.a[i] = 1/Branch.a[i]
	end
	Branch.fi[i] *= 3.14159265359/180
	Branch.smax[i] /= Sbase
end


if IData == 1
	@printf "------------------------------------------------------------------------------\n"
	@printf "---------------------------------INPUT DATA-----------------------------------\n"
	@printf "------------------------------------------------------------------------------\n"
	@printf "                                  BUS_DATA                                    \n"
	@printf "------------------------------------------------------------------------------\n"
	@printf "   Bus  Type    V0       Th0        Pd        Qd        Pg0        Qg0        \n" 
	@printf "------------------------------------------------------------------------------\n"
	for i in 1:nbus
		@printf "%5d" float(i)
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
	@printf " Branch From  To     r         x           z2         bshb        a        fi       smax \n" 
	@printf "-----------------------------------------------------------------------------------------\n"
	for i in 1:nbranch
		@printf "%5d" float(i)
		@printf "%5d" float(Branch.from[i])
	    @printf "%5d" float(Branch.to[i])
		@printf "%10.4f" float(Branch.r[i])
		@printf "%10.4f" float(Branch.x[i])
		@printf "%14.8f" float(Branch.z2[i])
		@printf "%10.4f" float(Branch.c[i])
		@printf "%10.4f" float(Branch.a[i])
		@printf "%10.4f" float(Branch.fi[i])
		@printf "%10.4f" float(Branch.smax[i])
		@printf "\n" 
	end
	@printf "-----------------------------------------------------------------------------------------\n"
	@printf "                                                                                         \n"
end


in_lines = [Int[] for i in 1:nbus] # indices of incoming branches
out_lines = [Int[] for i in 1:nbus] # indices of outgoing branches

#@show nbranch
#@show nbus

#a = zeros(Int, 1:nbranch)
#b = zeros(Int, 1:nbranch)
#@show(a)
#@show(b)
for i in 1:nbranch
#	a1 = int(Branch.from[i])
#	a = Int(Branch.from[i])

	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
#	@show(a)
    #@show(b)
    push!(out_lines[a],i)
    push!(in_lines[b],i)
#    @assert 1 <= Branch.to[i] <= nbus
end
#@show out_lines
#@show in_lines