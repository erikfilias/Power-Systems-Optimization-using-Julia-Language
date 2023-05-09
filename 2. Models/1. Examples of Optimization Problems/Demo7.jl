using JuMP
using Pajarito
using Pavito
using CPLEX
using Mosek
using Ipopt
baseprob = [0.3 0.4 0.3]
basedemand = [10.0 8.5 6.0]
basemu = [0.9 1.0 1.1]
baseϕ = [360.0 390.0 420.0]
prob = zeros(Float64, 729)
demand = zeros(Float64, 729)
mu=zeros(Float64, 4, 729)
ϕ =zeros(Float64, 729)
for i1 in 1:3
for i2 in 1:3
for i3 in 1:3
for i4 in 1:3
for i5 in 1:3
for i6 in 1:3
temp_scenario = 243*(i1-1) + 81*(i2-1) + 27*(i3-1) + 9*(i4-1) + 3*(i5-1) + i6
demand[temp_scenario] = basedemand[i1]
mu[1, temp_scenario] = basemu[i2]
mu[2, temp_scenario] = basemu[i3]
mu[3, temp_scenario] = basemu[i4]
mu[4, temp_scenario] = basemu[i5]
ϕ[temp_scenario] = baseϕ[i6]
prob[temp_scenario] = baseprob[i1] * baseprob[i2] * baseprob[i3] * baseprob[i4] * baseprob[i5] * baseprob[i6]
end
end
end
end
end
end

function generate_fullspace()
#set up model
# s1 = Model(solver=PajaritoSolver(rel_gap=0.0001, mip_solver=CplexSolver(), cont_solver=MosekSolver(MSK_IPAR_NUM_THREADS=0)))
s1 = Model(solver=PavitoSolver(rel_gap=0.0001, mip_solver=CplexSolver(), cont_solver=IpoptSolver()))

#sets for process
process = 1:4

#set for scenarios
scenarios = 1:729

#parameters
Qu = [2.3 2.8  2.0 3.2]
α = [80.0 100.0 70.0 110.0]
β =0.3* [90.0 80.0 100.0 72.0]
ψ = [1.0 1.0 1.0 1.0]
γ = 2.0*[50.0 50.0 50.0 50.0]


#first stage variables
@variable(s1, q[i in process]>=0.0)
@variable(s1, y[i in process], Bin)

#second stage variable
@variable(s1, R[i in process, s in scenarios]>=0.0)
@variable(s1, P[i in process, s in scenarios]>=0.0)
@variable(s1, z[i in process, s in scenarios], Bin)
@variable(s1, δ[s in scenarios]>=0.0 )

#original constraints in the first stage
@constraint(s1, e1[i in process], q[i]<= y[i] * Qu[i])

# Ax+g(y) ⩽0, g2(y)⩽0
@constraint(s1, c1[i in process, s in scenarios], z[i,s] <= y[i])
@constraint(s1, c2[i in process, s in scenarios], P[i,s]<= Qu[i] * z[i,s])
@constraint(s1, c3[s in scenarios], sum(P[i,s] for i in process) == demand[s] - δ[s])
@NLconstraint(s1, c4[i in process, s in scenarios], -log(1.0+R[i,s])* mu[i,s] + P[i,s] <= 0.0)
@constraint(s1, c5[i in process, s in scenarios], P[i,s] <= q[i])
#objective
@objective(s1, Min, sum(prob[s]*(ϕ[s]*δ[s] + sum(α[i]*y[i]+β[i]*q[i]+γ[i]*z[i,s] + ψ[i] * R[i,s] for i in process)) for s in scenarios  ) )
return s1
end

m = generate_fullspace()
print(m)
stat = solve(m)
println("objective value")
println(getobjectivevalue(m))
println("solver status")
println(stat)
println("First stage varialbes")
println("q:")
println(getvalue(getindex(m, :q)))
println("y:")
println(getvalue(getindex(m, :y)))
println("\nSecond stage varialbes")
println("R:")
println(getvalue(getindex(m, :R)))
println("P")
println(getvalue(getindex(m, :P)))
println("z")
println(getvalue(getindex(m, :z)))
println("δ")
println(getvalue(getindex(m, :δ)))`
