export AdditivePoissonDistortion

struct AdditivePoissonDistortion <: DistortionModel1D 
    rate::AbstractFloat 
    d::Poisson
end
function AdditivePoissonDistortion(λ::AbstractFloat)
    return AdditivePoissonDistortion(
        λ,
        Poisson(λ)
    )
end

function conditional_probabilities(model::AdditivePoissonDistortion, x::Integer, y::AbstractVector{<:Real})
    out = zeros(Float64, length(y))
    for i = 1:length(y)
        out[i] = (y[i] >= x) ? pdf(model.d, y[i]-x) : 0.0
    end
    return out 
end

function distortion(model::AdditivePoissonDistortion, x::Integer)::Number 
    return rand(model.d)
end