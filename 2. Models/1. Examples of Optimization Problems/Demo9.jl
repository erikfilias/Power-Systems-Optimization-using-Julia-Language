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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loading the necessary packages
#-------------------------------
using JuMP
using GLPKMathProgInterface
using GLPK
using LinearAlgebra
using Printf
using Random

# Master Problem Description
# --------------------------

#masterProblemModel = Model(solver=GLPKSolverMIP())
masterProblemModel = Model(with_optimizer(GLPK.Optimizer))

# Variable Definition
# ----------------------------------------------------------------
@variable(masterProblemModel, 0<= x[1:dimX]<= 1e6  , Int)
@variable(masterProblemModel, t<=1e6)

# Objective Setting
# -----------------
@objective(masterProblemModel, Max, t)
iC=1 # iC is the iteration counter

print(masterProblemModel)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Trying the entire benders decomposition in one step
while(true)
    global iC
    println("\n-----------------------")
    println("Iteration number = ", iC)
    println("-----------------------\n")
    println("The current master problem is")
    print(masterProblemModel)
    #statusMasterProblem = solve(masterProblemModel)
    optimize!(masterProblemModel)
    statusMasterProblem   = termination_status(masterProblemModel)
if statusMasterProblem == :Infeasible
    println("The problem is infeasible :-(")
    break
end

if statusMasterProblem == :Unbounded
    fmCurrent = M
    xCurrent=M*ones(dimX)
end


if statusMasterProblem == :Optimal
    fmCurrent = value.(t)
    xCurrent=Float64[]
    for i in 1:dimX
        push!(xCurrent,value.(x[i]))
    end
end

println("Status of the master problem is", statusMasterProblem,
        "\nwith fmCurrent = ", fmCurrent,
        "\nxCurrent = ", xCurrent)

# println("X1 \n")
subProblemModel = Model(solver=GLPKSolverLP())
# println("X2 \n")
cSub=b-A1*xCurrent
# println("X3 \n")
@variable(subProblemModel, u[1:dimU]>=0)
# println("X4 \n")

# @addConstraint(subProblemModel, constrRefSubProblem[j=1:size(A2,2)],sum(A2[i, j]*u[i] for i in [1:size(A2, 1)])>=c2[j])
@constraint(subProblemModel, constrRefSubProblem[j=1:size(A2,2),i=1:size(A2,1)],
sum(u[i]*A2[i,j])>=c2[j])
# println("X5 \n")
@show c1
@show xCurrent
Vec1 = 0
for i in 1:dimU
    Vec1 += cSub[i]*u[i]
end
# The second argument of @addConstraint macro, constrRefSubProblem[j=1:size(A2,2)] means that the j-th constraint is
# referenced by constrRefSubProblem[j].

@objective(subProblemModel, Min, dot(c1, xCurrent) + Vec1)
# @setObjective(subProblemModel, Min, dot(c1, xCurrent) + sum(cSub[i]*u[i] for i in [1:dimU]))

print("\nThe current subproblem model is \n", subProblemModel)

#statusSubProblem = solve(subProblemModel)
optimize!(subProblemModel)
statusSubProblem  = termination_status(subProblemModel)

fsxCurrent = objective_value(subProblemModel)

uCurrent = Float64[]

for i in 1:dimU
    push!(uCurrent, value.(u[i]))
end

γ=dot(b,uCurrent)

        println("Status of the subproblem is ", statusSubProblem,
        "\nwith fsxCurrent= ", fsxCurrent,
        "\nand fmCurrent= ", fmCurrent)

    if statusSubProblem == :Optimal &&  fsxCurrent == fmCurrent # we are done
        println("\n################################################")
        println("Optimal solution of the original problem found")
        println("The optimal objective value t is ", fmCurrent)
        println("The optimal x is ", xCurrent)
        println("The optimal v is ", getdual(constrRefSubProblem))
        println("################################################\n")
        break
    end

    if statusSubProblem == :Optimal && fsxCurrent < fmCurrent
        println("\nThere is a suboptimal vertex, add the corresponding constraint")
        cv= A1'*uCurrent - c1
        @constraint(masterProblemModel, t+sum(cv[i]*x[i] for i in 1:dimX) <= γ )
        println("t + ", cv, "ᵀ x <= ", γ)
    end

    if statusSubProblem == :Unbounded
        println("\nThere is an  extreme ray, adding the corresponding constraint")
        ce = A1'*uCurrent
        @constraint(masterProblemModel, sum(ce[i]*x[i] for i in 1:dimX) <= γ)
        println(ce, "ᵀ x <= ", γ)
    end

    iC=iC+1

end
