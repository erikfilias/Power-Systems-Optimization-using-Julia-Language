clearconsole()
using JuMP, GLPK, Printf
using Interact
using Plots

# Define some input data about the test system
# Maximum power output of generators
const GENERATION_MAX = [1000, 1000]
# Minimum power output of generators
const GENERATION_MIN = [0, 300]
# Incremental cost of generators
const COST_GENERATION = [50, 100]
# Incremental cost of wind generators
const COST_WIND = 50
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

    # Define decision variables
    # @variables(economic_dispatch, begin
    #     g[i=1:2]  # Thermal generation (MW).
    #     w >= 0  # Wind power (MW).
    # end)
    # g = Dict()
    # w = Dict()
    # for i in 1:2
    # 	@variable(economic_dispatch, g[i])
    # end
    # for i in 1:1
    # 	@variable(economic_dispatch, w[i] >= 0)
    # end
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


    # Return the optimal value of the objective function and its minimizers
    # as a NamedTuple.
    # gen = zeros(2)
    # wind_generation = zeros(2)
    # for i in 1:2
    #     gen[i] = value(g[i])
    #     wind_generation[i] = value(w)
    # end
    #
    # return (gen, wind_generation)
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


# for i in 1:2
println("   Generator[",1,"]: ", solution.generator1, " MW")
println("   Generator[2]: ", solution.generator2, " MW")
println("         Wind  : ", solution.wind_generation, " MW")
println("Wind spillage  : ", solution.wind_spillage, " MW")
println("----------------------------------")
println("Total cost: \$", solution.cost)


@manipulate for cost_of_wind in COST_WIND .* (1:0.1:3.5)
    solutions = Any[]
    cost_of_g1 = COST_GENERATION[1] .* (0.5:0.01:3.0)
    for c_g1 in cost_of_g1
        # update the incremental cost of the first generator at every iteration
        solution = solve_economic_dispatch(
            cost_of_thermal = [c_g1, COST_GENERATION[2]],
            cost_of_wind = cost_of_wind
        )
        push!(solutions, solution)
    end

    # Plot the outputs
    plot(
        # Plot the total cost
        plot(cost_of_g1, [sol.cost for sol in solutions],
            ylabel = "Total cost",
            ylims = (50000, 200000)
        ),
        # Plot the power output of Generator 1
        plot(cost_of_g1, [sol.generator1 for sol in solutions],
            ylabel = "Dispatch: G1",
            ylims = (0, 1100)
        ),
        # Plot the power output of Generator 2
        plot(cost_of_g1, [sol.generator2 for sol in solutions],
            ylabel = "Dispatch: G2",
            ylims = (0, 1600)
        ),
        # Plot the wind power output
        plot(cost_of_g1, [sol.wind_generation for sol in solutions],
            ylabel = "Dispatch: Wind",
            ylims = (0, 250)
        ),
        legend = false,
        xlabel = "Cost of G1"
    )
end
