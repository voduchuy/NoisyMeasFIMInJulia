# This subfolder is based on the two-state bursting gene transcription model. Adjust the following constants will change the parameters of the data-generating stochastic reaction network
const KON = 0.05 # Rate of gene activation
const KOFF = 0.015 # Rate of gene deactivation
const KR = 5.0 # Rate of mRNA sysnthesis
const GAMMA = 0.05 # Rate of mRNA degradation
const ALLOWED_MEASUREMENT_TIMES = Vector(0.0:1.0:400.0) # Vector of measurement times allowed

# All analyses in the `bursting_gene` subfolder are aimed at optimizing an experiment with the following characteristics:
# - Sampling times at `SAMPLING_TIME_COUNT` uniform time points.
# - Measuring `PER_TIMEPOINT_CELL_COUNT` i.i.d. single cells at each time point.
#The objective of the design is to maximize information by adjusting the timing between successive sampling times and the
#composition of measurement methods.
const PER_TIMEPOINT_CELL_COUNT = 1000
const SAMPLING_TIME_COUNT = 5

include(joinpath(@__DIR__, "./1d_distortions/Distortion1DModels.jl"))
include(joinpath(@__DIR__, "./utils/fim_utils.jl"))
include(joinpath(@__DIR__, "./utils/plot_utils.jl"))
using .Distortion1DModels
using .FimUtils
using .FIMPlotUtils

module ProbeBindingDefaultParameters
const PROBE_BINDING = 1.0E-1
const PROBE_UNBINDING = 1.0E-2
const CLUMPING_RATE = 1.0E-1
const CLUMP_DISSOLUTION = 5.0E-2
end