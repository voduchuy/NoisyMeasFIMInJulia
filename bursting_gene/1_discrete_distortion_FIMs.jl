include(joinpath(@__DIR__,"./common_settings.jl"))

using JLD2

# Truncation range for the domain and range state spaces of the probabilistic distortion operators
yrange = xrange = 0:400

kon, koff, α, γ = let f = load(joinpath(@__DIR__, "./results/bursting_parameters.jld2"))
    f["kon"], f["koff"], f["alpha"], f["gamma"]
end

distributions, sensitivities = let f = load(joinpath(@__DIR__, "./results/fsp_solutions.jld2"))
    f["rna_distributions"], f["rna_sensitivities"]
end

# Binomial distortion with uniform detection rate
binomial_uniform = BinomialDistortion(0.5) 
# Binomial distortion with state-dependent distortion rate, the more real RNA molecules, the more distorted
binomial_varying = BinomialDistortion() do x::Integer
    return 1.0 / (1.0 + 0.01 * x)
end
# Additive Poisson noise 
poisson_additive = AdditivePoissonDistortion(10.0)

#=
Compute the FIMs associated with noise-free mRNA counts 
=#
fim_exact = single_obs_fims(distributions, sensitivities)
jldsave(joinpath(@__DIR__,"./results/fim_exact.jld2"), fim=fim_exact)
#=
Compute FIMs associated single distorted observations based on the statistical models 
=#
single_distortion_models = Dict(["binomial"=>binomial_uniform, "binomial_state_dep"=>binomial_varying, "poisson"=>poisson_additive])
for (name, distortion_model) in single_distortion_models
    local C = dense_matrix(distortion_model, Vector(xrange), Vector(yrange))
    jldsave(joinpath(@__DIR__,"./results/distortion_matrix_$(name).jld2"), xrange=xrange, yrange=yrange, C=C)
    jldsave(joinpath(@__DIR__,"./results/fim_$(name).jld2"), fim=single_obs_fims(distributions, sensitivities, C))
end

#=
Compute the FIMs associated with composed distortions by chaining Binomial (uniform rate) -> Additive Poisson (in that order)
=# 
import LinearAlgebra: I as eye 
using LinearAlgebra: mul!
components = [poisson_additive, binomial_uniform]
C = eye(length(xrange))
for distortion in components 
    global C *= dense_matrix(distortion, Vector(xrange), Vector(yrange))
end
fim_composite = single_obs_fims(distributions, sensitivities, C);
jldsave(joinpath(@__DIR__,"./results/distortion_matrix_binomial_poisson.jld2"), xrange=xrange, yrange=yrange, C=C)
jldsave(joinpath(@__DIR__,"./results/fim_binomial_poisson.jld2"), fim=fim_composite)
