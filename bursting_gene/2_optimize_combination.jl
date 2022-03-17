include(joinpath(@__DIR__,"./common_settings.jl"))
using JLD2

const BUDGET_MAX = 1000.0
const MEASUREMENT_Δt = 30

#%%
struct MeasurementOption 
    title::String
    singlecell_fims::AbstractArray
    cost::Real 
end

function MeasurementOption(codename::String, title::String, cost::Real, θ::AbstractVector{<:Real})
    f = load(joinpath(@__DIR__,"results/fim_$(codename).jld2"))
    fims = f["fim"]    
    log_transform!(fims, θ)
    return MeasurementOption(
        title, 
        fims,
        cost 
    )
end

θ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    [f["kon"], f["koff"], f["alpha"], f["gamma"]]
end

#%%
method1 = MeasurementOption("exact","smFISH",1.0,θ);
method2 = MeasurementOption("binomial_state_dep","smFISH with random missing spots",0.5,θ);

#%%
obscount_max1 = Int(BUDGET_MAX/method1.cost)
obscount_max2 = Int(BUDGET_MAX/method2.cost)
objvalue_grid = zeros(obscount_max1+1,obscount_max2+1);
measurement_times = MEASUREMENT_Δt:MEASUREMENT_Δt:5*MEASUREMENT_Δt
for n1 in 0:obscount_max1
    for n2 in 0:obscount_max2
        if n1*method1.cost + n2*method2.cost > BUDGET_MAX
            continue  
        end
        𝔉1 = iidcombined_fim(method1.singlecell_fims, measurement_times, [n1])
        𝔉2 = iidcombined_fim(method2.singlecell_fims, measurement_times, [n2])
        𝔉 = 𝔉1 + 𝔉2 
        objvalue_grid[n1+1,n2+1] = compute_fim_functional(𝔉, :d)[]
    end 
end

#%%
opt_mixture = Tuple(argmax(objvalue_grid))
fim_mix_opt = iidcombined_fim(method1.singlecell_fims, measurement_times, [opt_mixture[1]]) + iidcombined_fim(method2.singlecell_fims, measurement_times, [opt_mixture[2]])

println("""
      The optimal mixture is $(opt_mixture[1]-1) $(method1.title) and
      $(opt_mixture[2]-1) $(method2.title) with D-opt value $(maximum(objvalue_grid))}.
      """)

jldsave(joinpath(@__DIR__,"results/opt_mixture.jld2"),
        opt_mixture = opt_mixture,
        obj_values = objvalue_grid,
        fim_mix_opt = fim_mix_opt
        )






