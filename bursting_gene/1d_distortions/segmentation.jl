export SegmentationDistortion 

using DSP: conv 

struct SegmentationDistortion <: DistortionModel1D
    rate::Real 
end

function (model::SegmentationDistortion)(pvec::AbstractVector{<:Real})::AbstractVector{<:Real}
    return model.rate*conv(pvec,pvec) + (1-model.rate)*[pvec;zeros(length(pvec)-1)]
end

function (model::SegmentationDistortion)(pvec::AbstractVector{<:Real}, svec::AbstractVector{<:Real})::AbstractVector{<:Real}
    return 2.0*model.rate * conv(svec, pvec) + (1.0-model.rate)*[svec;zeros(length(pvec)-1)]    
end