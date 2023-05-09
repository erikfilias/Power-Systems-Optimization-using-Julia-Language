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
system_name = "IEEE30"

# Adquisition DATA
include("SMC_dat.jl")

# Variables
@variable(m, Vsqr[Bus.busnum] >= 0)
@variable(m, th[Bus.busnum])
@variable(m, Pg[Bus.busnum])
@variable(m, Qg[Bus.busnum])

for i in 1:nbus
    set_start_value(Vsqr[i], Bus.V0[i]^2)
	set_start_value(th[i], Bus.Th0[i]*3.14159/180)
	set_start_value(Pg[i], Bus.Pg0[i])
	set_start_value(Qg[i], Bus.Qg0[i])
end
@variable(m, P[Branch.branchnum])
@variable(m, Q[Branch.branchnum])
@variable(m, Isqr[Branch.branchnum] >= 0)

for i in 1:nbranch
	set_start_value(P[i], 0)
	set_start_value(Q[i], 0)
	set_start_value(Isqr[i], 0)
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
	@NLconstraint(m, Pg[k] - Bus.Pd[k] + Vsqr[k]*Bus.gshb[k]
	 +	sum(P[i] for i in in_lines[k])
	 -	sum(P[i] + Branch.r[i]*Isqr[i] for i in out_lines[k]) == 0)
	#Q_balance_rule
	@NLconstraint(m, Qg[k] - Bus.Qd[k] + Vsqr[k]*Bus.bshb[k]
	 +	sum(Q[i] + Branch.c[i]*Vsqr[k] for i in in_lines[k])
	 -	sum(Q[i] - Branch.c[i]*Vsqr[k] + Branch.x[i]*Isqr[i]  for i in out_lines[k]) == 0)
end



for i in 1:nbranch
	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
	#Voltage difference
	@NLconstraint(m, Vsqr[a]*Branch.a[i]^2 - 2*(Branch.r[i]*P[i] + Branch.x[i]*Q[i]) - Branch.z2[i]*Isqr[i] - Vsqr[b] == 0)
	#Angle difference
	@NLconstraint(m, Branch.a[i]*sqrt(Vsqr[a])*sqrt(Vsqr[b])*sin(th[a] - th[b] + Branch.fi[i]) == Branch.x[i]*P[i] - Branch.r[i]*Q[i])
	#Current Magnitude
	@NLconstraint(m, Vsqr[b]*Isqr[i] == P[i]^2 + Q[i]^2)
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
		@constraint(m, Vsqr[i] == Bus.V0[i]^2)
	end
	if Bus.bustype[i] == 3
		@constraint(m, th[i] == Bus.Th0[i]*3.14159/180)
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
@printf "                                      BUS__RESULTS                                           \n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf "   Bus  Type    V0       Th0       Pg        Qg        Pd        Qd        gshb      bshb    \n"
@printf "---------------------------------------------------------------------------------------------\n"
for i in 1:nbus
    @printf "%5d" float(i)
    @printf "%5d" float(Bus.bustype[i])
    @printf "%10.4f" float(sqrt(value(Vsqr[i])))
    @printf "%10.4f" float(value(th[i])*180/3.14159265359)
    @printf "%10.4f" float(value(Pg[i])*Sbase)
    @printf "%10.4f" float(value(Qg[i])*Sbase)
    @printf "%10.4f" float(Bus.Pd[i]*Sbase)
    @printf "%10.4f" float(Bus.Qd[i]*Sbase)
    @printf "%10.4f" float(Sbase*Bus.gshb[i]*(value(Vsqr[i])))
    @printf "%10.4f" float(Sbase*Bus.bshb[i]*(value(Vsqr[i])))
    @printf "\n"
end
@printf "---------------------------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%35.4f" float(Sbase*sum(value(Pg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(value(Qg[i]) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Pd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.Qd[i] for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.gshb[i]*(value(Vsqr[i])) for i in 1:nbus))
@printf "%10.4f" float(Sbase*sum(Bus.bshb[i]*(value(Vsqr[i])) for i in 1:nbus))
@printf "\n"
@printf "---------------------------------------------------------------------------------------------\n"

@printf "---------------------------------------------------------------------------------------------\n"
@printf "                                     BRANCH__RESULTS                                         \n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf " Branch From  To  Pij[MW]   Pji[MW]   Qij[MW]   Qji[MW]    Pls[MW]   Qls[MVAr]   SLoad \n"
@printf "---------------------------------------------------------------------------------------------\n"

for i in 1:nbranch
	@printf "%5d" float(i)
    @printf "%5d" float(Branch.from[i])
    @printf "%5d" float(Branch.to[i])
    @printf "%10.4f" float(Sbase*(value(P[i])+Branch.r[i]*value(Isqr[i])))
    @printf "%10.4f" float(-Sbase*value(P[i]))
    @printf "%10.4f" float(Sbase*(value(Q[i])+Branch.x[i]*value(Isqr[i])-Branch.c[i]*value(Vsqr[Branch.from[i]])))
    @printf "%10.4f" float(-Sbase*(value(Q[i])+Branch.c[i]*value(Vsqr[Branch.to[i]])))
    @printf "%10.4f" float(Sbase*Branch.r[i]*value(Isqr[i]))
    @printf "%10.4f" float(Sbase*(Branch.x[i]*value(Isqr[i])-Branch.c[i]*value(Vsqr[Branch.from[i]])-Branch.c[i]*value(Vsqr[Branch.to[i]])))
    @printf "%10.2f" float(max((Sbase*(value(P[i])+Branch.r[i]*value(Isqr[i])))^2+(Sbase*(value(Q[i])+Branch.x[i]*value(Isqr[i])-Branch.c[i]*value(Vsqr[Branch.from[i]])))^2,(-Sbase*value(P[i]))^2+(-Sbase*(value(Q[i])+Branch.c[i]*value(Vsqr[Branch.to[i]])))^2)/(Branch.smax[i]*Sbase)^2*100)
    @printf "\n"
end

@printf "---------------------------------------------------------------------------------------------\n"
@printf "TOTAL"

@printf "%60.4f" float(Sbase*sum(Branch.r[i]*value(Isqr[i]) for i in 1:nbranch))
@printf "%10.4f" float(Sbase*sum(Branch.x[i]*value(Isqr[i]) - Branch.c[i]*value(Vsqr[Branch.from[i]]) - Branch.c[i]*value(Vsqr[Branch.to[i]]) for i in 1:nbranch))
@printf "\n"
@printf "---------------------------------------------------------------------------------------------\n"
@printf "\n"
