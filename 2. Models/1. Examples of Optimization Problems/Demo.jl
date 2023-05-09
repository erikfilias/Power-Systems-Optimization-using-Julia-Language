clearconsole()
#------------------------------------------Money Production Problem------------------------------------------
#----------How to have an amount of coins with a minimum total of weight
#-----------Value  Weight
# Pennies   1      2.5
# Nickels   5      5
# Dimes     10     2.268
# Quarters  25     5.670


# Packages(Librerias)
using JuMP
# using Cbc
using GLPK
using Printf

@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------GENESIS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"

# Model & Solver
# m = Model(with_optimizer(Cbc.Optimizer))
m = Model(with_optimizer(GLPK.Optimizer))

# Variables
@variable(m, pennies >= 0, Int)
@variable(m, nickels >= 0, Int)
@variable(m, dimes >= 0, Int)
@variable(m, quarters >= 0, Int)

# Constraints
@constraint(m, 1*pennies + 5*nickels + 10*dimes + 25*quarters == 99)

# Objective
@objective(m, Min, 2.5*pennies + 5*nickels + 2.268*dimes +5.670*quarters)

# Print Model
print(m)

# Initialization of the optimization
optimize!(m)
status = termination_status(m)
@printf "--------------------------------------------------------------------------------------\n"
@printf "-----------------------------------------RESULTS--------------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
println("Status of the Optimization: ", status)
println("Minimum weight: ", objective_value(m), "grams")
println("using:")
println(round(value.(pennies)), " pennies") #"round to cast as integer"
println(round(value.(nickels)), " nickels")
println(round(value.(dimes)), " dimes")
println(round(value.(quarters)), " quarters")
