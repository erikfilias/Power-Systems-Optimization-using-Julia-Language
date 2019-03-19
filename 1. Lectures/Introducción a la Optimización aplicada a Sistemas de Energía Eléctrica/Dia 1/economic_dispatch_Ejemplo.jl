clearconsole()
using JuMP, GLPK,Printf

# maxima generacion
const GENERATION_MAX = [1000, 1000]
# minima generacion
const GENERATION_MIN = [0, 300]
# Costos de generacion termicos
const COST_GENERATION = [50, 100]
# Costos de gen eolica
const COST_WIND = 10
# Demanda
const DEMAND = 1500
# WIND MAX
const WIND_MAX = 200;

function solve_economic_dispatch(;
    cost_of_thermal = COST_GENERATION,
    cost_of_wind = COST_WIND)

    economic_dispatch = Model(with_optimizer(GLPK.Optimizer))
    g = @variable(economic_dispatch, [1:2])
    @variable(economic_dispatch, w >= 0)

    @objective(economic_dispatch, Min, sum(cost_of_thermal[i]*g[i] for i in 1:2) + cost_of_wind*w)

    for i in 1:2
        @constraint(economic_dispatch, g[i] <= GENERATION_MAX[i])
        @constraint(economic_dispatch, g[i] >= GENERATION_MIN[i])
    end
    @constraints(economic_dispatch, begin
    # Max GEN WIND
     w <= WIND_MAX
    # Balance de energia
     sum(g[i] for i in 1:2) + w == DEMAND
    end)

    JuMP.optimize!(economic_dispatch)

    return (
        generator1 = value(g[1]),
        generator2 = value(g[2]),
        wind_generation = value(w),
        wind_spillage = WIND_MAX - value(w),
        cost = JuMP.objective_value(economic_dispatch),
        status = termination_status(economic_dispatch)
    )
end

solution = solve_economic_dispatch()
println("Status: ", solution.status)
println("Generator1: ", solution.generator1)
println("Generator2: ", solution.generator2)
println("Wind Gen: ", solution.wind_generation)
println("Wind Spillage: ", solution.wind_spillage)
println("Total cost: ", solution.cost)
