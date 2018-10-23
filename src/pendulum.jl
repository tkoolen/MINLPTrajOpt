module Pendulum

export
    PendulumParameters,
    MinEffort,
    MinTime,
    PendulumSwingUpProblem

using JuMP
using LinearAlgebra
using Parameters
using ..MINLPTrajOpt: sincosvar

@with_kw struct PendulumParameters
    m::Float64 = 1.0
    l::Float64 = 2.0
    g::Float64 = 9.81
    b::Float64 = 0.0#0.1
end

abstract type ObjectiveType end

struct MinEffort <: ObjectiveType end
struct MinTime <: ObjectiveType end

struct PendulumSwingUpProblem
    parameters::PendulumParameters
    model::JuMP.Model
    Δt::Vector{JuMP.Variable}
    sΔθ::Vector{JuMP.Variable}
    cΔθ::Vector{JuMP.Variable}
    Δθ::Vector{JuMP.Variable}
    sθ::Vector{JuMP.Variable}
    cθ::Vector{JuMP.Variable}
    θd::Vector{JuMP.Variable}
    θdd::Vector{JuMP.Variable}
    τ::Vector{JuMP.Variable}

    function PendulumSwingUpProblem(parameters::PendulumParameters, x0, solver;
            τmax::Number, # maximum torque
            N::Integer, # number of integration steps
            Δtmin::Number, # min time step
            Δtmax::Number, # max time step
            objectivetype::ObjectiveType,
            T::Union{Number, Nothing} = objectivetype isa MinTime ? nothing : N * Δtmax, # final time
            Δθmax::Number = 0.5, # max angle change during one integration time step (to control badness of small-angle approximation)
        )
        model = Model(solver=solver)

        # Time step settings
        @assert Δtmax >= Δtmin
        fixedstep = Δtmin == Δtmax
        if fixedstep
            T == N * Δtmax || throw(ArgumentError())
        end

        # Initial state
        θ0, θd0 = x0
        sθprev, cθprev = sincos(θ0)
        θdprev = θd0

        # Variable / expression vectors
        Δt = Vector{JuMP.Variable}(undef, N)
        sΔθ = Vector{JuMP.Variable}(undef, N)
        cΔθ = Vector{JuMP.Variable}(undef, N)
        Δθ = Vector{JuMP.Variable}(undef, N)
        sθ = Vector{JuMP.Variable}(undef, N)
        cθ = Vector{JuMP.Variable}(undef, N)
        θd = Vector{JuMP.Variable}(undef, N)
        θdd = Vector{JuMP.Variable}(undef, N)
        τ = Vector{JuMP.Variable}(undef, N)

        for i = 1 : N
            # Time steps
            Δt[i] = @variable model basename="Δt_{$i}"
            if fixedstep
                JuMP.fix(Δt[i], Δtmin)
                Δti = Δtmin
            else
                JuMP.setlowerbound(Δt[i], Δtmin)
                JuMP.setupperbound(Δt[i], Δtmax)
                Δti = Δt[i]
            end

            # Kinematics delta
            sΔθ[i], cΔθ[i] = sincosvar(model, "Δθ_{$i}", θmax=Δθmax, normconstraint=false)
            Δθ[i] = sΔθ[i] # first-order approximation

            # Absolute kinematics (from rotation matrix product)
            sθ[i], cθ[i] = sincosvar(model, "θ_{$i}", normconstraint=true)
            @NLconstraint model sθ[i] == sθprev * cΔθ[i] + cθprev * sΔθ[i]
            @NLconstraint model cθ[i] == cθprev * cΔθ[i] - sθprev * sΔθ[i]
            sθprev, cθprev = sθ[i], cθ[i]

            # Velocity, acceleration
            θd[i] = @variable model basename="θd_{$i}"
            @constraint model Δti * θd[i] == Δθ[i]
            θdd[i] = @variable model basename="θdd_{$i}"
            @constraint model Δti * θdd[i] == θd[i] - θdprev
            θdprev = θd[i]

            # Torque
            τ[i] = @variable model basename="τ_{$i}"
            setlowerbound(τ[i], -τmax)
            setupperbound(τ[i],  τmax)

            # Dynamics
            M = parameters.m * parameters.l^2
            mgl = parameters.m * parameters.g * parameters.l
            c = @NLexpression model parameters.b * θd[i] + mgl * sθ[i]
            @NLconstraint model M * θdd[i] + c == τ[i]
        end

        # Total time constraint
        if !fixedstep && T !== nothing
            @constraint model sum(Δt) == T
        end

        # Final state constraint
        θf = π
        θdf = 0.0
        @constraint model sθ[N] == sin(θf)
        @constraint model cθ[N] == cos(θf)
        @constraint model θd[N] == θdf

        # Objective
        if objectivetype isa MinEffort
            @objective model Min τ ⋅ τ
        elseif objectivetype isa MinTime
            @objective model Min sum(Δt)
        else
            throw(ArgumentError("Objective type not recognized"))
        end

        new(parameters, model, Δt, sΔθ, cΔθ, Δθ, sθ, cθ, θd, θdd, τ)
    end
end

end
