clearconsole()
using JuMP, GLPK, Printf

# Define some input data about the test system
# Maximum power output of generators
const GENERATION_MAX = [1000, 1000]
# Minimum power output of generators
const GENERATION_MIN = [0, 300]
# Incremental cost of generators
const COST_GENERATION = [50, 100]
# Incremental cost of wind generators
const COST_WIND = 10
# Total demand
const DEMAND = 1500
# Wind forecast
const WIND_MAX = 200;

"""
    solve_economic_dispatch(; cost_of_thermal::Vector, cost_of_wind)

Formulate and solve the economic dispatch problem given the cost of generation
for the two thermal generators and the cost of wind generation.
"""
function solve_economic_dispatch(;
        cost_of_thermal = COST_GENERATION,
        cost_of_wind = COST_WIND)

    economic_dispatch = Model(with_optimizer(GLPK.Optimizer))
    g = @variable(economic_dispatch, [1:2])
    @variable(economic_dispatch, w >= 0)

    # Define the objective function
    @objective(economic_dispatch, Min,
        sum(cost_of_thermal[i] * g[i] for i in 1:2) + cost_of_wind * w
    )

    # Define the constraint on the maximum and minimum power output of each generator.
    for i in 1:2
        @constraint(economic_dispatch, g[i] <= GENERATION_MAX[i])
        @constraint(economic_dispatch, g[i] >= GENERATION_MIN[i])
    end

    @constraints(economic_dispatch, begin
        # Define the constraint on the wind power injection
        w <= WIND_MAX
        # Define the power balance constraint
        sum(g[i] for i in 1:2) + w == DEMAND
    end)

    # Solve statement
    JuMP.optimize!(economic_dispatch)

    return (
        generator1 = value(g[1]),
        generator2 = value(g[2]),
        wind_generation = JuMP.value(w),
        wind_spillage = WIND_MAX - JuMP.value(w),
        cost = JuMP.objective_value(economic_dispatch),
        status = termination_status(economic_dispatch)
    )
end

# Solve the economic dispatch problem
solution = solve_economic_dispatch()
# solve_economic_dispatch()
@printf "--------------------------------------------------------------------------------------\n"
@printf "------------------------------------Results--Dispatch---------------------------------\n"
@printf "--------------------------------------------------------------------------------------\n"
println("Status of the Optimization: ", solution.status)
println("   Generator[",1,"]: ", solution.generator1, " MW")
println("   Generator[",2,"]: ", solution.generator2, " MW")
println("         Wind  : ", solution.wind_generation, " MW")
println("Wind spillage  : ", solution.wind_spillage, " MW")
println("----------------------------------")
println("Total cost: \$", solution.cost)
