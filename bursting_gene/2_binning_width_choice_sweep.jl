include(joinpath(@__DIR__,"./common_settings.jl"))
using JLD2
using LinearAlgebra: det

θ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    [f["kon"], f["koff"], f["alpha"], f["gamma"]]
end

distributions, sensitivities, t_meas = let f = load(joinpath(@__DIR__, "./results/fsp_solutions.jld2"))
    f["rna_distributions"], f["rna_sensitivities"], f["t_meas"]
end

#%%
bin_widths = 1:60
Δt_max = Int(t_meas[end]/SAMPLING_TIME_COUNT)

fim_dets = Vector{Vector{Float64}}()
optimal_dts = Vector{Int64}()
optimal_fims = Vector{Matrix}()
for (i,width) in enumerate(bin_widths)
    dmodel = Binning(width)
    C = dense_matrix(dmodel, Vector(0:400), Vector(0:400))
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

jldsave(joinpath(@__DIR__,"./results/binning_sweep.jld2"), 
    fim_dets=fim_dets,
    optimal_fims = optimal_fims,
    optimal_dts = optimal_dts,
    bin_widths = bin_widths,
    dt_max = Δt_max 
    )


