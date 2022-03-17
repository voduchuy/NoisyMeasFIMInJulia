module FIMPlotUtils
export confidence_ellipse!

using Plots 
using LinearAlgebra: eigen 

function confidence_ellipse!(p, fim::Matrix, num_sigma::Integer, parameter_ids::NTuple{2,<:Integer}, theta::AbstractVector{<:Real};label::String="",color, ls, lw::Real=2.0)    
    covmat = inv(fim)
    
    eigvals, eigvecs = let F = eigen(covmat[[parameter_ids[1], parameter_ids[2]], [parameter_ids[1], parameter_ids[2]]])
        F.values, F.vectors 
    end
    
    # Sort eigenvalues by descending order and re-order eigvecs correspondingly
    ids = sortperm(eigvals, rev=true)
    eigvals = eigvals[ids]
    eigvecs = eigvecs[:, ids]
    
    
    μ1 = theta[parameter_ids[1]]
    μ2 = theta[parameter_ids[2]]    
    σ1 = sqrt(eigvals[1])
    σ2 = sqrt(eigvals[2])
    a = num_sigma * σ1
    b = num_sigma * σ2    

    φ = atan(eigvecs[2, 1], eigvecs[1, 1])
    if φ < 0
        φ += 2 * π
    end
    
    rotation_matrix = [
        [cos(φ) -sin(φ)];
        [sin(φ) cos(φ)]
    ]    

    phi_grid = range(0, 2 * π, length=100)
    ellipse_x_r = a * cos.(phi_grid)
    ellipse_y_r = b * sin.(phi_grid)    

    r_ellipse = [ rotation_matrix*[x, y] for (x,y) in zip(ellipse_x_r,ellipse_y_r)]

    plot!(p, [x[1] for x in r_ellipse] .+ μ1, [x[2] for x in r_ellipse] .+ μ2, color=color, label=label, ls=ls)

    nothing 
end
end