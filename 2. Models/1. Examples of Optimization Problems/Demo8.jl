# Data for the problem
                   #---------------------
c1=[-1;-4]
c2=[-2; -3]
dimX=length(c1)
dimU=length(c2)
b=[-2;-3]
A1=[
    1 -3;
   -1 -3
   ]
A2=[
    1 -2;
   -1 -1
   ]
M=1000


# Loading the necessary packages
                    #-------------------------------
using JuMP
using LinearAlgebra
#using Compat.LinearAlgebra
using Printf
using Random
using CPLEX
# using Gurobi

                    # Master Problem Description
                    # --------------------------

# Model name

#masterProblemModel = Model(solver = CplexSolver())
masterProblemModel = Model(with_optimizer(CPLEX.Optimizer))
# masterProblemModel = Model(solver = GurobiSolver(Heuristics=0.0, Cuts = 0)) # If we want to add Benders lazy constraints
# in Gurobi, then we have to turn of Gurobi's own Cuts and Heuristics in the master problem

# Variable Definition (Only CplexSolver() works properly for these)
# ----------------------------------------------------------------
@variable(masterProblemModel,  x[1:dimX] >=0 , Int)
@variable(masterProblemModel, t)

# ***************ALTERNATIVE VARIABLE DEFINITION FOR GUROBI************
#If we replace the two lines above with the follwoing:
#@defVar(masterProblemModel,  0<= x[1:dimX] <= 1e6 , Int)
#@defVar(masterProblemModel, t <= 1e6)
# then all the solvers give the expected solution
#**********************************************************************

# Objective Setting
# -----------------
@objective(masterProblemModel, Max, t)

print(masterProblemModel)

stringOfBenderCuts=String[] # this is an array of strings which will contain all the
# Benders cuts added to be displayed later

# There are no constraints when we start, so we will add all the constraints in the
# form of Benders cuts lazily



function addBendersCut(cb)
    #***************************************************************************
    # First we store the master problem solution in conventional data structures
    println("----------------------------")
    println("ITERATION NUMBER = ", length(stringOfBenderCuts)+1)
    println("---------------------------\n")

    fmCurrent = value.(t)
    xCurrent=Float64[]
    for i in 1:dimX
        push!(xCurrent,value.(x[i]))
    end

    # Display the current solution of the master problem
    println("MASTERPROBLEM INFORMATION")
    println("-------------------------")
    println("The master problem that was solved was:")
    print(masterProblemModel)
    println("with ", length(stringOfBenderCuts), " added lazy constraints")
    println(stringOfBenderCuts)
    println("Current Value of x is: ", xCurrent)
    println("Current objective value of master problem, fmCurrent is: ", fmCurrent)
    println("\n")

    #************************************************************************

    # ========================================================================
    #                         Now we solve the subproblem
    println("X1 \n")
    #subProblemModel=Model(solver=CplexSolver())
    subProblemModel = Model(with_optimizer(CPLEX.Optimizer))
    println("X2 \n")
    #  subProblemModel = Model(solver=GurobiSolver(Presolve=0.0))

    cSub=b-A1*xCurrent
    println("X3 \n")
    @variable(subProblemModel, u[1:dimU]>=0)
    println("X4 \n")

    # @addConstraint(subProblemModel, constrRefSubProblem[j=1:size(A2,2)], sum(u[i]*A2[i,j] for i in [1:size(A2,1)])>=c2[j])
    @constraint(subProblemModel, constrRefSubProblem[j=1:size(A2,2),i=1:size(A2,1)], sum(u[i]*A2[i,j])>=c2[j])
    # for j in 1:size(A2,2)
    #     for i in 1:size(A2,1)
    #         @addConstraint(subProblemModel, sum(u[i]*A2[i,j])>=c2[j])
    #     end
    # end
    println("X5 \n")
    @show c1
    @show xCurrent
    Vec1 = 0
    for i in 1:dimU
        Vec1 += cSub[i]*u[i]
    end
    @objective(subProblemModel, Min, dot(c1, xCurrent) + Vec1)
    # @objective(subProblemModel, Min, dot(c1, xCurrent) + sum(cSub[i]*u[i] for i in [1:dimU]))
    println("X6 \n")
    println("The subproblem is being solved")

    #statusSubProblem = solve(subProblemModel)
    optimize!(subProblemModel)

    # We store the results achieved from the subproblem in conventional data structures

    fsxCurrent = objective_value(subProblemModel)

    uCurrent = Float64[]
    for i in 1:dimU
        push!(uCurrent, value.(u[i]))
    end

    # Display the solution corresponding to the subproblem

    println("SUBPROBLEM INFORMATION")
    println("----------------------")
    println("The subproblem that was solved was: ")
    print(subProblemModel)
    println("Current status of the subproblem is ", statusSubProblem)
    println("Current Value of u is: ", uCurrent) # JuMP will return an extreme ray
    # automatically (if the solver supports it), so we do not need to change the syntax
    println("Current Value of fs(xCurrent) is: ", fsxCurrent)
    println("\n")

    # ==========================================================================



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Now we check the status of the algorithm and add Benders cut when necessary
    γ=dot(b,uCurrent)



    if statusSubProblem == :Optimal &&  fsxCurrent==fmCurrent # we are done
        println("OPTIMAL SOLUTION OF THE ORIGINAL PROBLEM FOUND :-)")
        println("The optimal objective value t is ", fmCurrent)
        println("The optimal x is ", xCurrent)
        println("The optimal v is ", getdual(constrRefSubProblem))
        println("\n")
        return
    end

    println("-------------------ADDING LAZY CONSTRAINT----------------")
        if statusSubProblem == :Optimal && fsxCurrent < fmCurrent
        println("\nThere is a suboptimal vertex, add the corresponding constraint")
        cv= A1'*uCurrent - c1
        @lazyconstraint(cb, t+sum(cv[i]*x[i] for i in 1:dimX) <= γ )
        println("t + ", cv, "ᵀ x <= ", γ)
        push!(stringOfBenderCuts, string("t+", cv, "'x <=", γ))
    end

    if statusSubProblem == :Unbounded
        println("\nThere is an  extreme ray, adding the corresponding constraint")
        ce = A1'*uCurrent
        @lazyconstraint(cb, sum(ce[i]*x[i] for i in 1:dimX) <= γ)
        println(ce, "x <= ", γ)
        push!(stringOfBenderCuts, string(ce, "ᵀ x <= ", γ))
    end
    println("\n")
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end


addlazycallback(masterProblemModel, addBendersCut) # Telling the solver to use the
# callback function

#solve(masterProblemModel)
optimize!(masterProblemModel)
