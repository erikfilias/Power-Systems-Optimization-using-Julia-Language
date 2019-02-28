clearconsole()
# Packages(Librerias)
using JuMP, GLPK, Printf

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
# Model & Solver
m = Model(with_optimizer(GLPK.Optimizer))
# Variables
@variable(m, 0 <= x <= 10)
@variable(m, 0 <= y <=  5)
# Objective
@objective(m, Max, 10x + 8y )
# Constraints
@constraint(m,  4.5x + 1.5y <= 30 )
@constraint(m,  6x + 3y <= 48)
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
