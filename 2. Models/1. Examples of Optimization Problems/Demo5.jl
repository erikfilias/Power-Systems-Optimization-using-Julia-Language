using JuMP, AmplNLWriter
using CPLEX
using MosekTools
using GLPK
# using SCS

#MODEL CONSTRUCTION
#------------------

# sfLpModel = Model(solver = GLPKSolverLP())
#sfLpModel = Model(solver = CplexSolver())
# JuMP 0.19
# sfLpModel = Model(with_optimizer(CPLEX.Optimizer))
# sfLpModel = Model(with_optimizer(Mosek.Optimizer))
#sfLpModel = Model(with_optimizer(SCS.Optimizer))
 sfLpModel = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/Computer/ampl/bonmin"))

#INPUT DATA
#----------

c = [1; 3; 5; 2]

A= [
     1 1 9 5;
     3 5 0 8;
     2 0 6 13
    ]

b = [7; 3; 5]

m, n = size(A) # m = number of rows of A, n = number of columns of A

#VARIABLES
#---------

@variable(sfLpModel, x[1:n] >= 0) # Models x >=0

#CONSTRAINTS
#-----------

for i=1:m # for all rows do the following
    @constraint(sfLpModel, sum(A[i,j] * x[j] for j = 1:n) == b[i]) # the ith row
    # of A*x is equal to the ith component of b
end # end of the for loop

#OBJECTIVE
#---------
#typeof(objective_function(sfLpModel, QuadExpr))

@NLobjective(sfLpModel, Min, sum(c[j] * x[j]^2 for j = 1:n)) # minimize c'x
#sense = MOI.MIN_SENSE
#@objective(sfLpModel, sense, sum(c[j] * x[j]^2 for j = 1:n))

#THE MODEL IN A HUMAN-READABLE FORMAT
#------------------------------------
#@printf "-------------------------------------------------------------------------------------------------\n"
println("The optimization problem to be solved is:")
print(sfLpModel) # Shows the model constructed in a human-readable form

#@time begin
#status = solve(sfLpModel) # solves the model
optimize!(sfLpModel)
#end
#SOLVE IT AND DISPLAY THE RESULTS
#--------------------------------

#println("Objective value: ", getobjectivevalue(sfLpModel)) # getObjectiveValue(model_name) gives the optimum objective value
println("Objective value: ", objective_value(sfLpModel)) # getObjectiveValue(model_name) gives the optimum objective value

#println("Optimal solution is x = \n", getvalue(x)) # getValue(decision_variable) will give the optimum value
println("Optimal solution is x = \n", value.(x)) # getValue(decision_variable) will give the optimum value
                                                   # of the associated decision variable
