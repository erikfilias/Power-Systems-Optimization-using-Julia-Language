clearconsole()
using JuMP, Ipopt, Printf

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
# Model & Solver
m = Model(with_optimizer(Ipopt.Optimizer))

# Sistema a Simular
system_name = "IEEE14"

# Adquisition DATA
include("SMC_dat.jl")

# Variables
@variable(m, V[Bus.busnum] >= 0)
for i in 1:nbus
    if Bus.bustype[i] == 0
    	setvalue(V[i], Vnom)
    end
end
@variable(m, th[Bus.busnum])
for i in 1:nbus
    if Bus.bustype[i] != 3
    	setvalue(th[i], 0)
    end
end
@variable(m, P[Branch.branchnum])
for i in 1:nbranch
	setvalue(P[i], 0)
end
@variable(m, Q[Branch.branchnum])
for i in 1:nbranch
	setvalue(Q[i], 0)
end
@variable(m, Pg[Bus.busnum])
#for i in 1:nbus
#    if Bus.bustype[i] == 3
#    	setvalue(Pg[i], 0)
#    else
#    	setvalue(Pg[i], Bus.Pg0[i])
#    end
#end
@variable(m, Qg[Bus.busnum])

#@show Vsqr
#@show th
#@show Isqr
#@show P
#@show Q
#@show Pg
#@show Qg
@printf "-----------------------------------------------------------------------------------------\n"
@printf "                                                                                         \n"
@objective(m, Min, sum(Branch.r[i]*Isqr[i] for i=1:nbranch))
#end
for k in 1:nbus
	#P_balance_rule
	@constraint(m, Pg[k] - Bus.Pd[k] + Vsqr[k]*Bus.gshb[k] +
	#	sum{P[j,i], j in out_lines[k], i in in_lines[k]} -
		sum(P[i] for i in in_lines[k]) -
	#	sum{P[i,j] + Bus.r[i]*Isqr[i,j], i in out_lines[k], j in out_lines[k]} == 0 )
		sum(P[i] + Branch.r[i] * Isqr[i] for i in out_lines[k]) == 0)
	#Q_balance_rule
	@constraint(m, Qg[k] - Bus.Qd[k] + Vsqr[k]*Bus.bshb[k] +
	#	sum{P[j,i], j in out_lines[k], i in in_lines[k]} -
		sum(Q[i] + Branch.c[i] * Vsqr[k] for i in in_lines[k]) -
	#	sum{P[i,j] + Bus.r[i]*Isqr[i,j], i in out_lines[k], j in out_lines[k]} == 0 )
		sum(Q[i] - Branch.c[i] * Vsqr[k] + Branch.x[i] * Isqr[i] for i in out_lines[k]) == 0)
end

for i in 1:nbranch
	#Voltage_rule
	a = convert(Int, Branch.from[i])
	b = convert(Int, Branch.to[i])
	@constraint(m, Vsqr[a]*Branch.a[i]^2 -
		Vsqr[b]  == 2*(Branch.r[i]*P[i] + Branch.x[i]*Q[i]) + Branch.z2[i]*Isqr[i])
	#Angle_rule
	@NLconstraint(m, squareroot(Vsqr[a])*Branch.a[i]*squareroot(Vsqr[b])*sin(th[a] -
		th[b] + Branch.fi[i]) == Branch.x[i]*P[i] - Branch.r[i]*Q[i])
	#Current_rule
	@NLconstraint(m, Vsqr[b]*Isqr[i] == P[i]^2 +Q[i]^2)
end

#FIXED VARiABLE
for i in 1:nbus
	if Bus.bustype[i] == 3
		@constraint(m, th[i] == Bus.Th0[i]*3.14159265359/180)
	end
	if Bus.bustype[i] != 0
		@constraint(m, Vsqr[i] == Bus.V0[i]^2)
	end
	if Bus.bustype[i] != 3
		@constraint(m, Pg[i] == Bus.Pg0[i])
	end
	if Bus.bustype[i] == 0
		@constraint(m, Qg[i] == Bus.Qg0[i])
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
println("Utilidad: ", JuMP.objective_value(m))
println("Mesas  = ", JuMP.value(x))
println("Sillas = ", JuMP.value(y))
