include(joinpath(@__DIR__,"./common_settings.jl"))
using JLD2
using DataStructures
using Distributions: Categorical

mcsample_count = 10000

kon, koff, α, γ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    f["kon"], f["koff"], f["alpha"], f["gamma"]
end

distributions, sensitivities = let f = load(joinpath(@__DIR__, "./results/fsp_solutions.jld2"))
    f["rna_distributions"], f["rna_sensitivities"]
end

distortion_model = IntegratedIntensity(25.0, sqrt(25.0), 200.0, 400.0)
xrange = Vector(0:400)
yrange = Vector(range(distortion_model.μ_BG - 3*distortion_model.σ_BG,distortion_model.μ_probe*400+3*distortion_model.σ_BG,1000))
C = dense_matrix(distortion_model, xrange, yrange)
jldsave(joinpath(@__DIR__,"./results/distortion_matrix_flowcyt.jld2"), C=C, xrange=xrange, yrange=yrange)