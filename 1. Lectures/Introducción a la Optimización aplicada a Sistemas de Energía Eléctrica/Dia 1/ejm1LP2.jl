clearconsole()

using JuMP, GLPK

m = Model(with_optimizer(GLPK.Optimizer))

#c= [1 ; 2; 5]
c= [1  2 5]'
A= [-1 1 3 ; 1 3 -7]
b= [-5; 10]

@variable(m, x[1:3] >= 0)

@objective(m, Max, sum( c[i]* x[i] for i in 1:3))

@constraint(m, constraint[j in 1:2],  sum(A[j,i]*x[i] for i in 1:3) <= b[j] ) 
@constraint(m, bound, x[1] <= 10)

print(m)

#JuMP.optimize!(m)
optimize!(m)
status = termination_status(m)

println("Status of the Optimization: ", status)
println("Objective value: ", objective_value(m))
for i  in 1:3
 println("x[$i] =", value.(x[i]))
end
for j  in 1:2
 println("dual[$j] =", dual.(constraint[j]))
end
