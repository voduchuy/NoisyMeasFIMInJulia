export IntegratedIntensity

Base.@kwdef struct IntegratedIntensity <: DistortionModel1D
    μ_probe::AbstractFloat 
    σ_probe::AbstractFloat 
    μ_BG::AbstractFloat
    σ_BG::AbstractFloat 
end

function conditional_probabilities(model::IntegratedIntensity, x::Integer, y::AbstractVector{<:Real})::Vector{<:Real}
    μ_probe = model.μ_probe
    σ_probe = model.σ_probe 
    μ_BG = model.μ_BG 
    σ_BG = model.σ_BG 
    out = zeros(Float64, length(y))
    for i in 1:length(y)
        out[i] = exp(-(y[i] - x*μ_probe - μ_BG)^2.0/(2.0*(x*σ_probe^2.0 + σ_BG^2.0)))/sqrt(2*π*(x*σ_probe^2.0 + σ_BG^2.0))        
    end
    out
end

function distortion(model::IntegratedIntensity, x::Integer)
    return rand(Normal(model.μ_probe*x + model.μ_BG, sqrt(model.σ_probe^2*x + model.σ_BG^2)))
end