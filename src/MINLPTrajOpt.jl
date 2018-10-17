module MINLPTrajOpt

using JuMP
using LinearAlgebra
using Parameters

export
    PendulumParameters,
    PendulumSwingUpProblem

include("util.jl")
include("pendulum.jl")

end # module
