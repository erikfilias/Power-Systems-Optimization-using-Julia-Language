clearconsole()

using JuMP, GLPK

m = Model(with_optimizer(GLPK.Optimizer))

@variable(m,  0 <=x1 <= 10)
@variable(m,  x2 >= 0)
@variable(m,  x3 >= 0)

@objective(m, Max, x1 + 2x2 + 5x3)

@constraint(m,c1, -x1 + x2 + 3x3 <= -5)
@constraint(m,c2, x1 + 3x2 - 7x3 <= 10)

print(m)

#JuMP.optimize!(m)
optimize!(m)
status = termination_status(m)

println("Status of the Optimization: ", status)
println("Objective value: ", objective_value(m))
println("x1 = ", value.(x1))
println("x2 = ", value.(x2))
println("x3 = ", value.(x3))
