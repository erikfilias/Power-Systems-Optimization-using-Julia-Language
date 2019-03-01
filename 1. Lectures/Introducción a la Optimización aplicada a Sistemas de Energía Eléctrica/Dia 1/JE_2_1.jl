clearconsole()
# Packages(Librerias)
using JuMP, GLPK, Printf

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
# Model & Solver
m = Model(with_optimizer(GLPK.Optimizer))
# Variables
@variable(m, 0 <= x )
@variable(m, 0 <= y )
# Objective
@objective(m, Min, 11x + 9y )
# Constraints
@constraint(m, 1000000 <= (0.40)x + (0.32)y )
@constraint(m,  400000 <= (0.20)x + (0.40)y )
@constraint(m,  250000 <= (0.35)x + (0.20)y )
# Print Model
print(m)

# Initialization of the optimization
JuMP.optimize!(m)
status = termination_status(m)
@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------RESULTS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
println("Status of the Optimization: ", status)
println("Costo: ", JuMP.objective_value(m))
println("Petroleo crudo ligero = ", JuMP.value(x))
println("Petroleo crudo pesado = ", JuMP.value(y))
