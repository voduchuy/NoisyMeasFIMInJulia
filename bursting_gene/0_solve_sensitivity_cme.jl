using NumCME
using SparseArrays, LinearAlgebra
import DifferentialEquations as DE
import Sundials
using JLD2
include(joinpath(@__DIR__,"./common_settings.jl"))


function telegraph_model(kon, koff, kr, gamma)
    S = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
    a1 = propensity() do x, p
        p[1] * x[1]
    end
    a2 = propensity() do x, p
        p[2] * x[2]
    end
    a3 = propensity() do x, p
        p[3] * x[2]
    end
    a4 = propensity() do x, p
        p[4] * x[3]
    end
    return CmeModelWithSensitivity(S, [a1, a2, a3, a4], [kon, koff, kr, gamma])
end

model = telegraph_model(KON, KOFF, KR, GAMMA)
init_cond = forwardsens_initial_condition([[1, 0, 0]], [1.0], [[0.0] for i in 1:4])
fspmethod = AdaptiveForwardSensFspSparse(
    space_adapter = ForwardSensRStepAdapter(10, 10, true),
    ode_method = Sundials.CVODE_BDF(linear_solver = :GMRES)
)
sol = solve(model,
    init_cond,
    (0.0, ALLOWED_MEASUREMENT_TIMES[end]),
    fspmethod;
    saveat = ALLOWED_MEASUREMENT_TIMES,
    fsptol = 1.0E-6,
    odeatol = 1.0E-14,
    odertol = 1.0E-6, verbose = true)

rna_distributions = []
rna_sensitivities = []
for it in 1:length(sol)
    prna = Array(sum(sol[it].p, [1, 2]))
    Srna = [Array(sum(sol[it].S[ip], [1, 2])) for ip in 1:4]
    push!(rna_distributions, prna)
    push!(rna_sensitivities, Srna)
end

kon, koff, alpha, gamma = get_parameters(model)
jldsave(joinpath(@__DIR__,"./results/bursting_parameters.jld2"); kon, koff, alpha, gamma)
jldsave(joinpath(@__DIR__,"./results/fsp_solutions.jld2"); rna_distributions, rna_sensitivities, ALLOWED_MEASUREMENT_TIMES)



