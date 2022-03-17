include(joinpath(@__DIR__,"./common_settings.jl"))
using JLD2

const CONCENTRATION_LEVELS = [0.1, 1.0, 5.0, 10.0]
distortion_matrices = Dict{Float64,Matrix}()
for probe_concentration in CONCENTRATION_LEVELS
    println(probe_concentration)
    println("======================")
    distortion = ProbeBindingDistortion(
        ProbeBindingDefaultParameters.PROBE_BINDING,
        ProbeBindingDefaultParameters.PROBE_UNBINDING,
        ProbeBindingDefaultParameters.CLUMPING_RATE,
        ProbeBindingDefaultParameters.CLUMP_DISSOLUTION,
        probe_concentration
    )
    x = 0:400
    y = 0:400
    C = dense_matrix(distortion, x, y)
    distortion_matrices[probe_concentration] = C
end

jldsave(joinpath(@__DIR__,"results/smfish_probe_distortions.jld2"), 
                pdos = distortion_matrices, levels=CONCENTRATION_LEVELS)

#%%
θ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    [f["kon"], f["koff"], f["alpha"], f["gamma"]]
end

distributions, sensitivities, t_meas = let f = load(joinpath(@__DIR__, "./results/fsp_solutions.jld2"))
    f["rna_distributions"], f["rna_sensitivities"], f["t_meas"]
end

#%%
using LinearAlgebra: det 
Δt_max = Int(t_meas[end]/SAMPLING_TIME_COUNT)

fim_dets = Vector{Vector{Float64}}()
optimal_dts = Vector{Int64}()
optimal_fims = Vector{Matrix}()
for (i,level) in enumerate(CONCENTRATION_LEVELS)    
    C = distortion_matrices[level]
    singlecell_fims = single_obs_fims(distributions, sensitivities, C)
    log_transform!(singlecell_fims, θ)
    push!(fim_dets, 
        [
            det(
                iidcombined_fim(singlecell_fims, dt:dt:SAMPLING_TIME_COUNT*dt, [PER_TIMEPOINT_CELL_COUNT])
            )
            for dt in 1:Δt_max 
        ]
    )
    push!(optimal_dts, argmax(fim_dets[end]))
    push!(optimal_fims, let dt = optimal_dts[end]
                iidcombined_fim(singlecell_fims, dt:dt:SAMPLING_TIME_COUNT*dt, [PER_TIMEPOINT_CELL_COUNT])
            end
    )
end

jldsave(joinpath(@__DIR__,"./results/smfish_probe_sweep.jld2"), 
    fim_dets=fim_dets,
    optimal_fims = optimal_fims,
    optimal_dts = optimal_dts,
    levels=CONCENTRATION_LEVELS,
    dt_max = Δt_max 
    )
