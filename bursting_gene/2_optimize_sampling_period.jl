include(joinpath(@__DIR__,"./common_settings.jl"))
using JLD2
using DataStructures
using LinearAlgebra: det



θ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    [f["kon"], f["koff"], f["alpha"], f["gamma"]]
end

t_meas = let f = load(joinpath(@__DIR__, "./results/fsp_solutions.jld2"))
    f["t_meas"]
end

# %%
singlecell_fims = Dict()

for s in [
    "exact",
    "binomial",
    "poisson",
    "binomial_poisson",
    "flowcyt",
    "binomial_state_dep",
]
    f = load(joinpath(@__DIR__,"results/fim_$s.jld2")) 
    singlecell_fims[s] = f["fim"]    
    log_transform!(singlecell_fims[s], θ)    
end

#%% D-optimal sampling periods for different types of measurements
multicell_fims = Dict()

for meas_type in keys(singlecell_fims)
    multicell_fims[meas_type] = PER_TIMEPOINT_CELL_COUNT * singlecell_fims[meas_type]
end


dt_min = 1
dt_max = UInt32(floor(t_meas[end] / SAMPLING_TIME_COUNT))
dt_array = Vector(dt_min:dt_max)

multipletime_fims = Dict()
multipletime_fim_dets = Dict()

for meas_type in keys(multicell_fims)
    combined_fim = zeros(4,4,length(dt_array))
    det_comb_fim = zeros(length(dt_array))

    for i in 1:length(dt_array) 
        Δt = dt_array[i]
        combined_fim[:, :, i] = sum([ (@view multicell_fims[meas_type][:,:,j+1]) for j in Δt:Δt:SAMPLING_TIME_COUNT*Δt])
        det_comb_fim[i] = compute_fim_functional(combined_fim[:,:,i], :d)[]        
    end

    multipletime_fims[meas_type] = combined_fim
    multipletime_fim_dets[meas_type] = det_comb_fim
end

#%%
opt_rates = Dict()
for meas in keys(multipletime_fim_dets)
    Δt_opt = dt_array[argmax(multipletime_fim_dets[meas])]
    opt_rates[meas] = Δt_opt 
    println(
        """
        Optimal sampling period for $meas is $(Δt_opt) min with D-opt=$(multipletime_fim_dets[meas][Δt_opt])."""
    )
end 

jldsave(joinpath(@__DIR__,"./results/opt_sampling_periods.jld2"), opt_rates=opt_rates, multipletime_fims=multipletime_fims, multipletime_fim_dets=multipletime_fim_dets)

