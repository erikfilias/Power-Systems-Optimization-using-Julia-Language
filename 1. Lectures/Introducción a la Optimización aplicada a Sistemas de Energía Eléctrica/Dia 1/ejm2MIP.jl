clearconsole()

using JuMP, GLPK,  Printf

m = Model(with_optimizer(GLPK.Optimizer))

@variable(m, 0<= x1 <= 10)
@variable(m,  x2 >= 0, Int)
@variable(m,  x3, Bin)

@objective(m, Max, x1 + 2x2 +5x3)

@constraint(m, -x1 + x2 + 3x3 <= -5)
@constraint(m, x1 + 3x2 - 7x3 <= 10)

optimize!(m)
print(m)
status = termination_status(m)
@printf " ------ Solucion -------\n"
println("Estado de la optimización:  ",status)
println("Función objetivo: ",objective_value(m))
println("x1 = ",value(x1))
println("x2 = ",value(x2))
println("x3 = ",value(x3))
