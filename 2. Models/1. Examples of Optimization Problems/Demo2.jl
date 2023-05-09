clearconsole()
# Packages(Librerias)
using JuMP
# using Clp
using GLPK
using Printf

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"

# Model & Solver
# m = Model(with_optimizer(Clp.Optimizer))
m = Model(with_optimizer(GLPK.Optimizer))

# Variables
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

# Constraints
@constraint(m, 1x + 5y <= 3.0 )

# Objective
@objective(m, Max, 5x + 3*y )

# Print Model
print(m)

# Initialization of the optimization
optimize!(m)
status = termination_status(m)
@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------RESULTS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
println("Status of the Optimization: ", status)
println("Objective value: ", objective_value(m))
println("x = ", value.(x))
println("y = ", value.(y))
