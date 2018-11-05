module SO3

export
    SO3Parameters,
    SO3Problem

using JuMP
using LinearAlgebra
using Parameters
using Rotations

function Rotations.Quat(model::JuMP.Model, basename::AbstractString; normconstraint=true, θmax=π)
    w = @variable model basename=basename*"w"
    x = @variable model basename=basename*"x"
    y = @variable model basename=basename*"y"
    z = @variable model basename=basename*"z"
    q = Quat(w, x, y, z, false)

    # Add bounds to help the solver out:
    setlowerbound.((w, x, y, z), -1.0)
    setupperbound.((w, x, y, z),  1.0)

    # Avoid unnecessary branches due to Quat double cover.
    # Helps a lot.
    # setlowerbound(w, 0.0)

    if normconstraint
        @NLconstraint model w^2 + x^2 + y^2 + z^2 == 1
    end

    if θmax < π
        @constraint model w >= cos(θmax / 2)
        for var in [x, y, z]
            @constraint model sin(-θmax / 2) <= var <= sin(θmax / 2)
        end
    end
    q
end

# type piracy:
JuMP.getvalue(q::Quat{<:JuMP.AbstractJuMPScalar}) = Quat(getvalue(q.w), getvalue(q.x), getvalue(q.y), getvalue(q.z), false)

function quatmul(model::JuMP.Model, q1::Quat, q2::Quat, q12::Quat)
    @NLconstraint model q12.w == q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
    @NLconstraint model q12.x == q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y
    @NLconstraint model q12.y == q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x
    @NLconstraint model q12.z == q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    model
end

function varvec3(model::JuMP.Model, basename::AbstractString, index::Integer)
    x = @variable model basename="$(basename)_{$index,x}"
    y = @variable model basename="$(basename)_{$index,y}"
    z = @variable model basename="$(basename)_{$index,z}"
    [x, y, z]
end

@with_kw struct SO3Parameters
    Ixx::Float64 = 1.0
    Iyy::Float64 = 1.0
    Izz::Float64 = 1.0
end

struct SO3Problem
    parameters::SO3Parameters
    model::JuMP.Model
    q0::Quat{Float64}
    Δq::Vector{Quat{JuMP.Variable}}
    q::Vector{Quat{JuMP.Variable}}
    ω::Vector{Vector{JuMP.Variable}}
    ωd::Vector{Vector{JuMP.Variable}}
    τ::Vector{Vector{JuMP.Variable}}

    function SO3Problem(parameters::SO3Parameters, q0::Quat, ω0::AbstractVector, solver;
            τmax::Number, # maximum torque
            N::Integer, # number of integration steps
            Δt::Number, # time step
            Δθmax::Number = 0.5, # max angle change during one integration time step (to control badness of small-angle approximation)
        )
        model = Model(solver=solver)

        # Initial state
        q0 = principal_value(q0)
        qprev = q0
        ωprev = ω0

        # Variable / expression vectors
        Δq = Vector{Quat{JuMP.Variable}}(undef, N)
        q = Vector{Quat{JuMP.Variable}}(undef, N)
        ω = Vector{Vector{JuMP.Variable}}(undef, N)
        ωd = Vector{Vector{JuMP.Variable}}(undef, N)
        τ = Vector{Vector{JuMP.Variable}}(undef, N)

        for i = 1 : N
            # Kinematics delta
            Δq[i] = Quat(model, "Δq_{$i}", θmax=Δθmax, normconstraint=false)

            # Absolute kinematics
            q[i] = Quat(model, "q_{$i}", normconstraint=true)
            quatmul(model, qprev, Δq[i], q[i])
            qprev = q[i]

            # Angular velocity
            ω[i] = varvec3(model, "ω", i)
            @constraint model 2 * Δt * ω[i] .== [Δq[i].x, Δq[i].y, Δq[i].z] # from small-angle approximation of quaternion -> axis-angle map

            # Angular acceleration
            ωd[i] = varvec3(model, "ωd", i)
            @constraint model Δt * ωd[i] .== ω[i] - ωprev
            ωprev = ω[i]

            # Torque
            τ[i] = varvec3(model, "τ", i)
            @constraint model τ[i] ⋅ τ[i] <= τmax^2

            # Dynamics (Euler's equations)
            Ixx = parameters.Ixx
            Iyy = parameters.Iyy
            Izz = parameters.Izz
            @NLconstraint model Ixx * ωd[i][1] + (Izz - Iyy) * ω[i][2] * ω[i][3] == τ[i][1]
            @NLconstraint model Iyy * ωd[i][2] + (Ixx - Izz) * ω[i][3] * ω[i][1] == τ[i][2]
            @NLconstraint model Izz * ωd[i][3] + (Iyy - Ixx) * ω[i][1] * ω[i][2] == τ[i][3]
        end

        # Final state constraint
        # @constraint model q[N].w == 1.0
        # @constraint model q[N].x == 0.0
        # @constraint model q[N].y == 0.0
        # @constraint model q[N].z == 0.0
        # @constraint model ω[N] .== 0.0

        # Objective
        qw = [q[i].w for i = 1 : N]
        qv = [[q[i].x, q[i].y, q[i].z] for i = 1 : N]
        wτ = 0.1
        γf = 1.0#0.2
        γ = γf^(1/(N - 1))
        @objective model Min sum(γ^(i - 1) * (wτ^2 * τ[i] ⋅ τ[i] + qv[i]⋅ qv[i]) for i = 1 : N) / N

        new(parameters, model, q0, Δq, q, ω, ωd, τ)
    end
end


end # module
