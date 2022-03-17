export BinomialDistortion

struct BinomialDistortion{RT} <: DistortionModel1D
    detection_rate::RT
end
detection_rate(model::BinomialDistortion{<:Any},x::Integer) = model.detection_rate(x)
detection_rate(model::BinomialDistortion{<:AbstractFloat},x::Integer) = model.detection_rate

function conditional_probabilities(model::BinomialDistortion, x::Integer, y::AbstractVector{<:Real})
    out = zeros(Float64, length(y))
    for i = 1:length(y)
        out[i] = pdf(Binomial(x, detection_rate(model,x)), y[i])
    end
    return out 
end

function distortion(model::DistortionModel1D, x::Integer)::Number 
    return rand(Binomial(x, detection_rate(model,x)))
end

