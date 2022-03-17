module FimUtils 
export single_obs_fims, log_transform!, iidcombined_fim, compute_fim_functional
using LinearAlgebra: det, eigmin

using EllipsisNotation

function compute_fim_functional(fims, func::Symbol)
    if func == :d 
        functional = det 
    elseif func == :e 
        functional = eigmin 
    else 
        throw(ArgumentError("Requested functional not implemented.")) 
    end

    out = zeros(Float64, size(fims)[3:ndims(fims)]...)
    for ind in CartesianIndices(out)
        out[ind] = functional(fims[:,:,ind])
    end
    out
end

function log_transform!(fims, parameters)
    parametercount = size(fims, 1)
    if size(fims, 2) ≠ parametercount
        throw(DimensionMismatch("User-supplied FIMs have invalid dimensions. First and second dimensions of the FIM array must match and are equal to number of parameters."))
    end
    for i in 1:parametercount
        for j in 1:parametercount
            fims[i,j,..] .*= parameters[i]*parameters[j]*log(10.0)*log(10.0)
        end
    end
end

function single_obs_fims(distributions, sensitivities, C=nothing)
    nt = length(distributions)
    parametercount = length(sensitivities[1])
    fims = zeros(Float64, parametercount, parametercount, nt)    
    for i in 1:nt
        n = length(distributions[i])
        if C === nothing 
            p = distributions[i]
            S = [sensitivities[i][ip] for ip in 1:parametercount]
        else
            p = (@view C[:, 1:n]) * distributions[i]
            S = [(@view C[:, 1:n]) * sensitivities[i][ip] for ip in 1:parametercount]
        end
        
        for ip in 1:parametercount
            for jp in 1:ip
                fims[ip, jp, i] = sum(S[ip] .* S[jp] ./ max.(p, 1.0E-14))
            end
        end
        for ip in 1:parametercount
            for jp in ip+1:parametercount
                fims[ip, jp, i] = fims[jp, ip, i]
            end
        end
    end
    return fims
end

"""
    iidcombined_fim(singlecell_fims, time_ids, sample_counts)

Compute the Fisher Information Matrix associated with a dataset that consists of i.i.d. observations at multiple timepoints given 
an array `singlecell_fims[1:np,1:np,1:nt]` of single-cell FIMs (last dimension is the time dimension), indices of observation times `time_ids`, and vector `sample_counts` of observation count per time point. 
"""
function iidcombined_fim(singlecell_fims, time_ids, sample_counts)
    if length(sample_counts) == 1
        sample_counts = repeat(Vector(sample_counts), length(time_ids))
    elseif length(sample_counts) ≠ length(time_ids)
        throw(DimensionMismatch("List of time indices and sample counts have different lengths."))
    end
    paramcount = size(singlecell_fims, 1)
    out = zeros(paramcount, paramcount)
    if isempty(time_ids)
        throw(ArgumentError("Empty list of measurement times is not acceptable."))
    end
    for it in 1:length(time_ids)
        out .+= sample_counts[it]*(@view singlecell_fims[:,:,time_ids[it]])
    end
    out 
end
end