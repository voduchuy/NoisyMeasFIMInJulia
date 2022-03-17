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

function simulate_observations(p, n, model::IntegratedIntensity)
    d = Categorical(p)
    xsamples = rand(d, n)
    ysamples = zeros(Float64, length(xsamples))
    @simd for i in 1:length(xsamples)
        ysamples[i] = distortion(model, xsamples[i])
    end
    return ysamples
end

fims = zeros(4,4,length(distributions))
mcsc = 1.0/mcsample_count
for i in 1:length(distributions)
    println(i)
    p = distributions[i]
    p[p.<0] .= 0.0
    p ./= sum(p)
    S = sensitivities[i]
    ysamples = simulate_observations(p, mcsample_count, distortion_model)    
    C = dense_matrix(distortion_model, Vector(0:400), ysamples)
    
    Cp = (@view C[:,1:length(p)])*p
    CS = [(@view C[:,1:length(p)])*s for s in S]
    for ip in 1:4 
        for jp in ip:4 
            fims[ip,jp,i] = mcsc*sum((CS[ip]./max.(Cp, 1.0E-14)).*(CS[jp]./max.(Cp, 1.0E-14)))
        end
        for jp in 1:ip-1
            fims[ip,jp,i] = fims[jp,ip,i]
        end
    end    
end
jldsave(joinpath(@__DIR__,"./results/fim_flowcyt.jld2"), fim=fims)