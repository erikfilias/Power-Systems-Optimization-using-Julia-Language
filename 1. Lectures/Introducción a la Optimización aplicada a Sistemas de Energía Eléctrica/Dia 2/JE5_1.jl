clearconsole()
using JuMP, Ipopt, Printf, LinearAlgebra

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
# Model & Solver
m = Model(with_optimizer(Ipopt.Optimizer))

Vnom= 1.00

# Sistema a Simular
system_name = "IEEE3"

# Adquisition DATA
include("SMC_dat.jl")

# Variables
@variable(m, V[Bus.busnum] >= 0)
for i in 1:nbus
    set_start_value(V[i], Bus.V0[i])
end
@variable(m, th[Bus.busnum])
for i in 1:nbus
    set_start_value(th[i], Bus.Th0[i]*3.14159/180)
end
@variable(m, Pde[Branch.branchnum])
for i in 1:nbranch
	set_start_value(Pde[i], 0)
end
@variable(m, Qde[Branch.branchnum])
for i in 1:nbranch
	set_start_value(Qde[i], 0)
end
@variable(m, Ppara[Branch.branchnum])
for i in 1:nbranch
	set_start_value(Pde[i], 0)
end
@variable(m, Qpara[Branch.branchnum])
for i in 1:nbranch
	set_start_value(Qde[i], 0)
end
@variable(m, Pg[Bus.busnum])
for i in 1:nbus
   	set_start_value(Pg[i], Bus.Pg0[i])
end
@variable(m, Qg[Bus.busnum])
for i in 1:nbus
   	set_start_value(Qg[i], Bus.Qg0[i])
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
#end
for k in 1:nbus
	#P_balance_rule
	@NLconstraint(m, Pg[k] - Bus.Pd[k] + V[k]^2*Bus.gshb[k]
	 -	sum(Ppara[i] for i in in_lines[k])
	 -	sum(Pde[i]  for i in out_lines[k]) == 0)
	#Q_balance_rule
	@NLconstraint(m, Qg[k] - Bus.Qd[k] + V[k]^2*Bus.bshb[k]
	 -	sum(Qpara[i] for i in in_lines[k])
	 -	sum(Qde[i] for i in out_lines[k]) == 0)
end



for i in 1:nbranch
	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
	#Pde
	@NLconstraint(m, Pde[i] == Branch.g[i]*Branch.a[i]^2*V[a]^2
	-Branch.a[i]*V[a]*V[b]*Branch.g[i]*cos(th[a]-th[b]+Branch.fi[i])
	-Branch.a[i]*V[a]*V[b]*Branch.b[i]*sin(th[a]-th[b]+Branch.fi[i]))
	#Qde
	@NLconstraint(m, Qde[i] == -(Branch.b[i]+Branch.c[i])*Branch.a[i]^2*V[a]^2
	-Branch.a[i]*V[a]*V[b]*Branch.g[i]*sin(th[a]-th[b]+Branch.fi[i])
	+Branch.a[i]*V[a]*V[b]*Branch.b[i]*cos(th[a]-th[b]+Branch.fi[i]))
	#Ppara
	@NLconstraint(m, Ppara[i] == Branch.g[i]*V[b]^2
	-Branch.a[i]*V[a]*V[b]*Branch.g[i]*cos(th[a]-th[b]+Branch.fi[i])
	+Branch.a[i]*V[a]*V[b]*Branch.b[i]*sin(th[a]-th[b]+Branch.fi[i]))
	#Qpara
	@NLconstraint(m, Qpara[i] == -(Branch.b[i]+Branch.c[i])*V[b]^2
	+Branch.a[i]*V[a]*V[b]*Branch.g[i]*sin(th[a]-th[b]+Branch.fi[i])
	+Branch.a[i]*V[a]*V[b]*Branch.b[i]*cos(th[a]-th[b]+Branch.fi[i]))
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
		@constraint(m, V[i] == Bus.Vgen[i])
	end
	if Bus.bustype[i] == 3
		@constraint(m, th[i] == Bus.Th0[i]*3.14159/180)
	end
end


# Print Model
print(m)

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
    @printf "%10.4f" float(value(V[i]))
    @printf "%10.4f" float(value(th[i])*180/3.14159)
    @printf "%10.4f" float(value(Pg[i])*Sbase)
    @printf "%10.4f" float(value(Qg[i])*Sbase)
    @printf "%10.4f" float(Bus.Pd[i]*Sbase)
    @printf "%10.4f" float(Bus.Qd[i]*Sbase)
    @printf "%10.4f" float(Sbase*Bus.gshb[i]*(value(V[i]))^2)
    @printf "%10.4f" float(Sbase*Bus.bshb[i]*(value(V[i]))^2)
    @printf "\n"
end
@printf "---------------------------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%35.4f" float(Sbase*sum(value(Pg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(value(Qg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Pd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Qd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.gshb[i]*(value(V[i]))^2 for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.bshb[i]*(value(V[i]))^2 for i in 1:nbus))
@printf "\n"
@printf "---------------------------------------------------------------------------------------------\n"

@printf "-----------------------------------------------------------------------------\n"
@printf "                               BRANCH__RESULTS                               \n"
@printf "-----------------------------------------------------------------------------\n"
@printf " Branch From  To  Pij[MW]   Pji[MW]   Qij[MW]   Qji[MW]    Pls[MW] Qls[MVAr] \n"
@printf "-----------------------------------------------------------------------------\n"

for i in 1:nbranch
    @printf "%5d" float(i)
    @printf "%5d" float(Branch.from[i])
    @printf "%5d" float(Branch.to[i])
    @printf "%10.4f" float(Sbase*value(Pde[i]))
    @printf "%10.4f" float(-Sbase*value(Ppara[i]))
    @printf "%10.4f" float(Sbase*value(Qde[i]))
    @printf "%10.4f" float(-Sbase*(value(Qpara[i])))
    @printf "%10.4f" float(Sbase*(value(Pde[i])+value(Ppara[i])))
    @printf "%10.4f" float(Sbase*(value(Qde[i])+value(Qpara[i])))
    @printf "\n"
end

@printf "-----------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%60.4f" float(Sbase*sum((value(Pde[i])+value(Ppara[i])) for i in 1:nbranch))
@printf "%10.4f" float(Sbase*sum((value(Qde[i])+value(Qpara[i])) for i in 1:nbranch))
@printf "\n"
