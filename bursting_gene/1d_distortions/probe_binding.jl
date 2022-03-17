using NumCME
using StaticArrays: @MVector
using Sundials: CVODE_BDF

export ProbeBindingDistortion

Base.@kwdef struct ProbeBindingDistortion <: DistortionModel1D
    probe_binding_rate::Real
    probe_unbinding_rate::Real
    clumping_rate::Real
    clump_dissolution::Real
    probe_concentration::Real
    texposure::Real
    cmemodel::CmeModel
end

function ProbeBindingDistortion(binding::Real, unbinding::Real, clumping::Real, clump_breakup::Real, probe_concentration::Real, texposure::Real=3E2)
    S = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
    θ = [binding, unbinding, clumping, clump_breakup, probe_concentration]
    λ1 = propensity() do x, p
        p[1] * p[5] * x[1]
    end
    λ2 = propensity() do x, p
        p[2] * x[2]
    end
    λ3 = propensity() do x, p
        p[3] * p[5]
    end
    λ4 = propensity() do x, p
        p[4] * x[3]
    end
    return ProbeBindingDistortion(
        probe_binding_rate=binding,
        probe_unbinding_rate=unbinding,
        clumping_rate=clumping,
        clump_dissolution=clump_breakup,
        probe_concentration=probe_concentration,
        texposure=texposure,
        cmemodel=CmeModel(S, [λ1, λ2, λ3, λ4], θ)
    )
end

function conditional_probabilities(model::ProbeBindingDistortion, x::Integer, y::AbstractVector{<:Real})
    println(x)

    p0 = FspVectorSparse([@MVector [x, 0, 0]], [1.0])
    fspmethod = AdaptiveFspSparse(
        ode_method=CVODE_BDF(linear_solver=:GMRES),
        space_adapter=RStepAdapter(10, 10, true)
    )
    sol = solve(
        model.cmemodel,
        p0,
        (0.0, model.texposure),
        fspmethod
        ;
        saveat=[model.texposure],
        fsptol=1.0E-6,
        verbose=false
    )
    pend = sol[end].p

    # The distributions of visualized RNA count and false spot count are independent
    p_visualized_rna = Array(sum(pend, [1, 3]))
    p_falsespots = Array(sum(pend, [1, 2]))
    pconv = conv(p_visualized_rna, p_falsespots)

    ydict = Dict([yi => i for (i, yi) in enumerate(y)])
    out = zeros(Float64, length(y))
    @simd for i in 1:length(pconv)
        id = get(ydict, i - 1, 0)
        if id ≠ 0
            @inbounds out[id] = pconv[i]
        end
    end
    return out
end