clearconsole()
using JuMP, Ipopt, Printf, LinearAlgebra

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
# Model & Solver
# solver = IpoptSolver()
# m = Model(solver = solver)
m = Model(with_optimizer(Ipopt.Optimizer))

Vnom= 1.00

# Sistema a Simular
system_name = "IEEE14"

# Adquisition DATA
include("SMC_dat.jl")

# Variables
@variable(m, e[Bus.busnum])
@variable(m, f[Bus.busnum])
@variable(m, Pg[Bus.busnum])
@variable(m, Qg[Bus.busnum])

for i in 1:nbus
    set_start_value(e[i], Vnom)
	set_start_value(f[i], 0)
	set_start_value(Pg[i], Bus.Pg0[i])
	set_start_value(Qg[i], Bus.Qg0[i])
end
@variable(m, Pde[Branch.branchnum])
@variable(m, Qde[Branch.branchnum])
@variable(m, Ppara[Branch.branchnum])
@variable(m, Qpara[Branch.branchnum])

for i in 1:nbranch
	set_start_value(Pde[i], 0)
	set_start_value(Qde[i], 0)
	set_start_value(Ppara[i], 0)
	set_start_value(Qpara[i], 0)
end

#@show Vsqr
#@show th
#@show Isqr
#@show P
#@show Q
#@show Pg
#@show Qg
@printf "-----------------------------------------------------------------------------------------\n"
@printf "                                                                                         \n"
@objective(m, Min, sum(Pg[i] for i=1:nbus if Bus.bustype[i] == 3))
# @objective(m, Min, sum(1e6*(Pde[i]+Ppara[i]) for i=1:nbranch))
#end
for k in 1:nbus
	#P_balance_rule
	@NLconstraint(m, Pg[k] - Bus.Pd[k]
	 -	sum(Ppara[i] for i in in_lines[k])
	 -	sum(Pde[i]  for i in out_lines[k]) == 0)
	#Q_balance_rule
	@NLconstraint(m, Qg[k] - Bus.Qd[k] + (e[k]^2+f[k]^2)*Bus.bshb[k]
	 -	sum(Qpara[i] for i in in_lines[k])
	 -	sum(Qde[i] for i in out_lines[k]) == 0)
end



for i in 1:nbranch
	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
	#Pde
	@NLconstraint(m, Pde[i] == Branch.g[i]*Branch.a[i]^2*(e[a]^2+f[a]^2)
	-Branch.a[i]*Branch.g[i]*(e[a]e[b]+f[a]f[b])
	+Branch.a[i]*Branch.b[i]*(e[a]f[b]-e[b]f[a]))
	#Qde
	@NLconstraint(m, Qde[i] == -(Branch.b[i]+Branch.c[i])*Branch.a[i]^2*(e[a]^2+f[a]^2)
	+Branch.a[i]*Branch.g[i]*(e[a]f[b]-e[b]f[a])
	+Branch.a[i]*Branch.b[i]*(e[a]e[b]+f[a]f[b]))
	#Ppara
	@NLconstraint(m, Ppara[i] == Branch.g[i]*(e[b]^2+f[b]^2)
	-Branch.a[i]*Branch.g[i]*(e[a]e[b]+f[a]f[b])
	-Branch.a[i]*Branch.b[i]*(e[a]f[b]-e[b]f[a]))
	#Qpara
	@NLconstraint(m, Qpara[i] == -(Branch.b[i]+Branch.c[i])*(e[b]^2+f[b]^2)
	-Branch.a[i]*Branch.g[i]*(e[a]f[b]-e[b]f[a])
	+Branch.a[i]*Branch.b[i]*(e[a]e[b]+f[a]f[b]))
end

#FIXED VARiABLE
for i in 1:nbus
	if Bus.bustype[i] != 3
		@constraint(m, Pg[i] == Bus.Pg0[i])
	end
	if Bus.bustype[i] == 0
		@constraint(m, Qg[i] == Bus.Qg0[i])
	end
	if Bus.bustype[i] != 0
		@constraint(m, (e[i]^2+f[i]^2) == Bus.V0[i]^2)
	end
	if Bus.bustype[i] == 3
		@constraint(m, f[i] == e[i]*tan(Bus.Th0[i]*3.14159/180))
	end
end


# Print Model
# print(m)

# Initialization of the optimization
JuMP.optimize!(m)
status = termination_status(m)
@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------RESULTS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
println("Status of the Optimization: ", status)

println("Objective value: ", JuMP.objective_value(m)*Sbase)
@printf "---------------------------------------------------------------------------------------------\n"
@printf "                                         BUS__RESULTS                                        \n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf "   Bus  Type    V0       Th0       Pg        Qg        Pd        Qd        gshb      bshb    \n"
@printf "---------------------------------------------------------------------------------------------\n"
for i in 1:nbus
    @printf "%5d" float(i)
    @printf "%5d" float(Bus.bustype[i])
    @printf "%10.4f" float(sqrt(value(e[i])^2+value(f[i])^2))
    @printf "%10.4f" float(atan(value(f[i])/value(e[i]))*180/3.14159265359)
    @printf "%10.4f" float(value(Pg[i])*Sbase)
    @printf "%10.4f" float(value(Qg[i])*Sbase)
    @printf "%10.4f" float(Bus.Pd[i]*Sbase)
    @printf "%10.4f" float(Bus.Qd[i]*Sbase)
    @printf "%10.4f" float(Sbase*Bus.gshb[i]*(value(e[i])^2+value(f[i])^2))
    @printf "%10.4f" float(Sbase*Bus.bshb[i]*(value(e[i])^2+value(f[i])^2))
    @printf "\n"
end
@printf "---------------------------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%35.4f" float(Sbase*sum(value(Pg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(value(Qg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Pd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Qd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.gshb[i]*(value(e[i])^2+value(f[i])^2) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.bshb[i]*(value(e[i])^2+value(f[i])^2) for i in 1:nbus))
@printf "\n"
@printf "---------------------------------------------------------------------------------------------\n"

@printf "-----------------------------------------------------------------------------\n"
@printf "                               BRANCH__RESULTS                               \n"
@printf "-----------------------------------------------------------------------------\n"
@printf " Branch From  To  Pij[MW]   Pji[MW]   Qij[MW]   Qji[MW]    Pls[MW] Qls[MVAr] \n"
@printf "-----------------------------------------------------------------------------\n"

Pde = zeros(nbranch)
Ppara = zeros(nbranch)
Qde = zeros(nbranch)
Qpara = zeros(nbranch)
for i in 1:nbranch
	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
	Pde[i] = (
	Branch.a[i]^2*Branch.g[i]*value(e[a]^2+f[a]^2)
	-Branch.a[i]*Branch.g[i]*value(e[a]e[b]+f[a]f[b])
	+Branch.a[i]*Branch.b[i]*value(e[a]f[b]-e[b]f[a])
	)

	Ppara[i] = (
	Branch.g[i]*(value(e[b]^2+f[b]^2))
	-Branch.a[i]*Branch.g[i]*(value(e[a]e[b]+f[a]f[b]))
	-Branch.a[i]*Branch.b[i]*(value(e[a]f[b]-e[b]f[a]))
	)

	Qde[i] = (
	-Branch.a[i]^2*(Branch.b[i]+Branch.c[i])*(value(e[a]^2+f[a]^2))
	+Branch.a[i]*Branch.g[i]*(value(e[a]f[b]-e[b]f[a]))
	+Branch.a[i]*Branch.b[i]*(value(e[a]e[b]+f[a]f[b]))
	)

	Qpara[i] = (
	-(Branch.b[i]+Branch.c[i])*(value(e[b]^2+f[b]^2))
	-Branch.a[i]*Branch.g[i]*(value(e[a]f[b]-e[b]f[a]))
	+Branch.a[i]*Branch.b[i]*(value(e[a]e[b]+f[a]f[b]))
	)
end

for i in 1:nbranch
    @printf "%5d" float(i)
    @printf "%5d" float(Branch.from[i])
    @printf "%5d" float(Branch.to[i])
    @printf "%10.4f" float(Sbase*Pde[i])
    @printf "%10.4f" float(Sbase*Ppara[i])
    @printf "%10.4f" float(Sbase*Qde[i])
    @printf "%10.4f" float(Sbase*Qpara[i])
    @printf "%10.4f" float(Sbase*(Pde[i]+Ppara[i]))
    @printf "%10.4f" float(Sbase*(Qde[i]+Qpara[i]))
    @printf "\n"
end

@printf "-----------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%60.4f" float(Sbase*sum(Pde[i]+Ppara[i] for i in 1:nbranch))
@printf "%10.4f" float(Sbase*sum(Qde[i]+Qpara[i] for i in 1:nbranch))
@printf "\n"
