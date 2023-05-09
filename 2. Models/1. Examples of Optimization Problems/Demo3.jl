clearconsole()
using JuMP, AmplNLWriter
#using Clp
# using GLPKMathProgInterface
using GLPK
using Ipopt
using CPLEX

mutable struct Edge
	from; to; cost; capacity;
end
edges = [Edge(1,2,1,0.5), Edge(1,3,2,0.4), Edge(1,4,3,0.6),
		 Edge(2,5,2,0.3), Edge(3,5,2,0.6), Edge(4,5,2,0.5)]
node = 1:5

#solver = ClpSolver()
#solver = GLPKSolverLP()
# solver = GLPKSolverMIP()
#solver = IpoptSolver()
# solver = CplexSolver()
#m = Model(solver = solver)
#m = Model()
m = Model(with_optimizer(GLPK.Optimizer))
# m = Model(with_optimizer(CPLEX.Optimizer))
# m = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/Computer/ampl/cplex"))
# m = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/Computer/ampl/gurobi"))

@variable(m, 0 <= flow[e in edges] <= e.capacity)
@constraint(m,  sum(flow[e] for e in edges if e.to == 5) == 1)
#@constraint(m, flowcon[n=2:4], sum(flow[e] for e in edges if e.to == node) == sum(flow[e] for e in edges if e.from == node) )

@objective(m, Min, sum(e.cost * flow[e] for e in edges) )

# print(m)
#status = solve(m)
#optimize!(m, with_optimizer(GLPK.Optimizer))
optimize!(m)
if termination_status(m) == MOI.OPTIMAL
	println("Objective value: ", objective_value(m))
	println("Flow = ", value.(flow))
elseif termination_status(m) == MOI.TIME_LIMIT && has_values(m)
	 suboptimal_solution = value.(flow)
	 suboptimal_objective = objective_value(m)
else
	  error("The model was not solved correctly.")
end
	#println("Objective value: ", getobjectivevalue(m))
	#println("Flow = ", getvalue(flow))
	#println("y = ", getvalue(y))
