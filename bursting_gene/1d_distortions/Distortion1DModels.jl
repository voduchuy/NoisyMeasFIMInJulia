module Distortion1DModels
using Distributions
using DSP: conv

export DistortionModel1D, dense_matrix, conditional_probabilities, distortion

abstract type DistortionModel1D end

function conditional_probabilities(model::DistortionModel1D, x::Integer, y::AbstractVector{<:Real})
    zeros(Float64, length(y))
end
function dense_matrix(model::DistortionModel1D, xrange::AbstractVector{<:Real}, yrange::AbstractVector{<:Real})::Matrix
    out = zeros(Float64, length(yrange), length(xrange))
    @simd for i in 1:length(xrange)
        x = xrange[i]
        out[:, i] .= conditional_probabilities(model, x, yrange)
    end
    out
end
function distortion(model::DistortionModel1D, x::Integer)::Number
    throw(ErrorException("Please implement the method for $(typeof(model))."))
end
function (model::DistortionModel1D)(pvec::AbstractVector{<:Real})::AbstractVector{<:Real}
    throw(ErrorException("Please implement the method for $(typeof(model))."))
end
function (model::DistortionModel1D)(pvec::AbstractVector{<:Real}, svec::AbstractVector{<:Real})::AbstractVector{<:Real}
    throw(ErrorException("Please implement the method for $(typeof(model))."))
end

include("binomial_distortions.jl")
include("poisson_additive.jl")
include("binning_uniform.jl")
include("integrated_intensity.jl")
include("segmentation.jl")
include("probe_binding.jl")
end