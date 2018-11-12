module MINLPTrajOpt

export
    TrajOptProblem,
    MinEffort,
    MinTime

using JuMP
using LinearAlgebra
using RigidBodyDynamics
using TrigonometricPolynomials
using RigidBodyDynamics.CustomCollections: SegmentedVector
using RigidBodyDynamics: upper, lower

const JointSegmentedVector{T} = SegmentedVector{JointID, T, Base.OneTo{JointID}, Vector{T}}
using TrigonometricPolynomials: Var

include("problem.jl")
include(joinpath("objectives", "min_effort.jl"))
include(joinpath("objectives", "min_time.jl"))

# TODO: remove:
include("util.jl")
include("pendulum.jl")
include("acrobot.jl")
include("so3.jl")

end # module

