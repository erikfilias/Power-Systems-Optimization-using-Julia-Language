clearconsole()
using JuMP, AmplNLWriter  # Need to say it whenever we use JuMP
using CPLEX
# using GLPKMathProgInterface # Loading the GLPK module for using its solver
using GLPK
# using Mosek

#MODEL CONSTRUCTION
#--------------------

#myModel = Model(solver=GLPKSolverLP())
# myModel = Model(with_optimizer(GLPK.Optimizer))
# myModel = Model(with_optimizer(CPLEX.Optimizer))
# myModel = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/Computer/ampl/cplex"))
myModel = Model(with_optimizer(AmplNLWriter.Optimizer, "/Users/Computer/ampl/xpress"))

# myModel = Model(solver = CplexSolver(CPXPARAM_Preprocessing_Relax=0))
# myModel = Model(solver = MosekSolver(Problemtype(4)))
# Name of the model object. All constraints and variables of an optimization problem are associated
# with a particular model object. The name of the model object does not have to be myModel, it can be yourModel too! The argument of Model,
# solver=GLPKsolverLP() means that to solve the optimization problem we will use GLPK solver.

#VARIABLES
#---------

# A variable is modelled using @defVar(name of the model object, variable name and bound, variable type)
# Bound can be lower bound, upper bound or both. If no variable type is defined, then it is treated as
#real. For binary variable write Bin and for integer use Int.

@variable(myModel, x >= 0) # Models x >=0

# Some possible variations:
# @defVar(myModel, x, Binary) # No bound on x present, but x is a binary variable now
# @defVar(myModel, x <= 10) # This one defines a variable with lower bound x <= 10
# @defVar(myModel, 0 <= x <= 10, Int) # This one has both lower and upper bound, and x is an integer

@variable(myModel, y >= 0) # Models y >= 0

#OBJECTIVE
#---------

@objective(myModel, Max, x + y) # Sets the objective to be minimized. For maximization use Max

#CONSTRAINTS
#-----------

@constraint(myModel, 3x + 2y <= 5) # Adds the constraint x + y <= 1

#THE MODEL IN A HUMAN-READABLE FORMAT
#------------------------------------
println("The optimization problem to be solved is:")
print(myModel) # Shows the model constructed in a human-readable form

#SOLVE IT AND DISPLAY THE RESULTS
#--------------------------------
@time begin
#status = solve(myModel) # solves the model
optimize!(myModel)
end

#println("Objective value: ", getobjectivevalue(myModel)) # getObjectiveValue(model_name) gives the optimum objective value
println("Objective value: ", objective_value(myModel)) # getObjectiveValue(model_name) gives the optimum objective value
#println("x = ", getvalue(x)) # getValue(decision_variable) will give the optimum value of the associated decision variable
println("x = ", value.(x)) # getValue(decision_variable) will give the optimum value of the associated decision variable
#println("y = ", getvalue(y))
println("y = ", value.(y))
