clearconsole()
using LinearAlgebra
# bus = readdlm("IEEE30a.bus")
#
# branch = readdlm("IEEE30a.branch")

R =[1 7 0; 0 0 2; 1 0 5]
RanR = rank(R)
@show R
@show RanR


n = 5
p = 4
m = 3
A=
[0.7511 -0.1357   0.7955  -0.4567 0.1356
-0.6670 -0.3326   0.1657  -0.5519 -0.9367
 1.5894 -0.1302  -0.4313  -0.4875  0.4179]

B=
[-0.09520 -0.28056 -1.33978 0.6506
 -0.8581  -0.3518   1.2788  1.5114
 -0.5925  1.3477    0.1589  0.03495]

c=[0.3468,0.8687,0.1200,0.5024,0.2884]

d=[0.2017,0.2712,0.4997,0.9238]

f = [0.1716,0.3610,0.0705]

using JuMP, AmplNLWriter
using GLPK
#using CPLEX

#sfMipModel = Model(solver = CplexSolver())
# sfMipModel = Model(with_optimizer(CPLEX.Optimizer))
sfMipModel = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/jchavez/ampl/xpress"))
# sfMipModel = Model(solver = GLPKSolverMIP())

@variable(sfMipModel, x[1:n] >=0)
@variable(sfMipModel, y[1:p] >= 0, Int)

@objective(sfMipModel, Min, sum(c[i] * x[i] for i in 1:n) + sum(d[i] * y[i] for i in 1:p))

for i in [1:m]
    @constraint(sfMipModel, sum(A[i,j] * x[j] for j in [1:n]) + sum(B[i,j] * y[j] for j in [1:p]) .== f[i])
end

print(sfMipModel, "\n")
#statusMipModel = solve(sfMipModel)
optimize!(sfMipModel)
print("Status of the problem is ", termination_status(sfMipModel), "\n")

#if statusMipModel == :Optimal
if termination_status(sfMipModel) == MOI.OPTIMAL
    print("Optimal objective value = ", objective_value(sfMipModel), "\nOptimal x = ", value.(x), "\nOptimal y = ", value.(y))
end
print("Optimal objective value = ", objective_value(sfMipModel), "\nOptimal x = ", value.(x), "\nOptimal y = ", value.(y))
