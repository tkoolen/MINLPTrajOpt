module Acrobot

export
    AcrobotParameters,
    MinEffort,
    MinTime,
    AcrobotSwingUpProblem

using JuMP
using LinearAlgebra
using Parameters
using ..MINLPTrajOpt: sincosvar

@with_kw struct AcrobotParameters:
    m1::Float64 = 1.0
    l1::Float64 = 1.0
    m2::Float64 = 1.0
    l2::Float64 = 1.0
    g::Float64 = 9.81
    I1::Float64 = m1 * l1^2
    I2::Float64 = m2 * l2^2 
end

abstract type ObjectiveType end

struct MinEffort <: ObjectiveType end
struct MinTime <: ObjectiveType end

struct AcrobotSwingUpProblem
    parameters::AcrobotParameters
    model::JuMP.Model
    Δt::Vector{JuMP.Variable}
    sΔθ1::Vector{JuMP.Variable}
    cΔθ1::Vector{JuMP.Variable}
    Δθ1::Vector{JuMP.Variable}
    sθ1::Vector{JuMP.Variable}
    cθ1::Vector{JuMP.Variable}
    θ1d::Vector{JuMP.Variable}
    θ1dd::Vector{JuMP.Variable}
    sΔθ2::Vector{JuMP.Variable}
    cΔθ2::Vector{JuMP.Variable}
    Δθ2::Vector{JuMP.Variable}
    sθ2::Vector{JuMP.Variable}
    cθ2::Vector{JuMP.Variable}
    θ2d::Vector{JuMP.Variable}
    θ2dd::Vector{JuMP.Variable} 
    τ::Vector{JuMP.Variable}

    function AcrobotSwingUpProblem(parameters::AcrobotParameters, x0, solver;
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
        θ10, θ1d0, θ20, θ2d0 = x0
        sθ1prev, cθ1prev = sincos(θ10)
        sθ2prev, cθ2prev = sincos(θ20)
        θ1dprev = θ1d0
        θ2dprev = θ2d0

        # Variable / expression vectors
        Δt = Vector{JuMP.Variable}(undef, N)
        sΔθ1 = Vector{JuMP.Variable}(undef, N)
        cΔθ1 = Vector{JuMP.Variable}(undef, N)
        Δθ1 = Vector{JuMP.Variable}(undef, N)
        sθ1 = Vector{JuMP.Variable}(undef, N)
        cθ1 = Vector{JuMP.Variable}(undef, N)
        θ1d = Vector{JuMP.Variable}(undef, N)
        θ1dd = Vector{JuMP.Variable}(undef, N)
        sΔθ2 = Vector{JuMP.Variable}(undef, N)
        cΔθ2 = Vector{JuMP.Variable}(undef, N)
        Δθ2 = Vector{JuMP.Variable}(undef, N)
        sθ2 = Vector{JuMP.Variable}(undef, N)
        cθ2 = Vector{JuMP.Variable}(undef, N)
        θ2d = Vector{JuMP.Variable}(undef, N)
        θ2dd = Vector{JuMP.Variable}(undef, N)
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
            sΔθ1[i], cΔθ1[i] = sincosvar(model, "Δθ1_{$i}", θmax=Δθmax, normconstraint=false)
            sΔθ2[i], cΔθ2[i] = sincosvar(model, "Δθ2_{$i}", θmax=Δθmax, normconstraint=false)
            Δθ1[i] = sΔθ1[i] # first-order approximation
            Δθ2[i] = sΔθ2[i] # first-order approximation

            # Absolute kinematics (from rotation matrix product)
            sθ1[i], cθ1[i] = sincosvar(model, "θ1_{$i}", normconstraint=true)
            @NLconstraint model sθ1[i] == sθ1prev * cΔθ1[i] + cθ1prev * sΔθ1[i]
            @NLconstraint model cθ1[i] == cθ1prev * cΔθ1[i] - sθ1prev * sΔθ1[i]
            sθ1prev, cθ1prev = sθ1[i], cθ1[i]
            sθ2[i], cθ2[i] = sincosvar(model, "θ2_{$i}", normconstraint=true)
            @NLconstraint model sθ2[i] == sθ2prev * cΔθ2[i] + cθ2prev * sΔθ2[i]
            @NLconstraint model cθ2[i] == cθ2prev * cΔθ2[i] - sθ2prev * sΔθ2[i]
            sθ2prev, cθ2prev = sθ2[i], cθ2[i]

            # Velocity, acceleration
            θ1d[i] = @variable model basename="θ1d_{$i}"
            @constraint model Δti * θ1d[i] == Δθ1[i]
            θ1dd[i] = @variable model basename="θ1dd_{$i}"
            @constraint model Δti * θ1dd[i] == θ1d[i] - θ1dprev
            θ1dprev = θ1d[i]
            θ2d[i] = @variable model basename="θ2d_{$i}"
            @constraint model Δti * θ2d[i] == Δθ2[i]
            θ2dd[i] = @variable model basename="θ2dd_{$i}"
            @constraint model Δti * θ2dd[i] == θ2d[i] - θ2dprev
            θ2dprev = θ2d[i]

            # Torque
            τ[i] = @variable model basename="τ_{$i}"
            setlowerbound(τ[i], -τmax)
            setupperbound(τ[i],  τmax)

            # Dynamics
            c1 = parameters.I1 + parameters.I2 + parameters.m2 * parameters.l1^2
            c2 = parameters.m2 * parameters.l1 * parameters.l2
            c3 = parameters.g * parameters.l1 * (parameters.m1 + parameters.m2)
            c4 = parameters.g * parameters.m2 * parameters.l2
            n1 = @NLexpression model (c1 + 2 * c2 * cθ2[i]) * θ1dd[i]
            n2 = @NLexpression model (parameters.I2 + c2 * cθ2[i]) * θ2dd[i]
            n3 = @NLexpression model 2 * c2 * sθ2[i] * θ2d[i] * θ1d[i]
            n4 = @NLexpression model c2 * sθ2[i] * θ2d[i]^2
            n5 = @expression model c3 * sθ1[i]
            n6 = @expression model c4 * (sθ1[i]*cθ2[i]+cθ1[i]*sθ2[i])
            n7 = @NLexpression model (parameters.I2 + c2 * cθ2[i]) * θ1dd[i]
            n8 = @expression model parameters.I2 * θ2dd[i]
            n9 = @NLexpression model c2 * sθ2[i] * θ1d[i]^2

            @NLconstraint model n1 + n2 - n3 - n4 + n5 + n6 == 0
            @NLconstraint model n7 + n8 + n9 + n6 == τ[i]
        end

        # Total time constraint
        if !fixedstep && T !== nothing
            @constraint model sum(Δt) == T
        end

        # Final state constraint (swing up)
        θ1f = π
        θ1df = 0.0
        θ2f = 0.0
        θ2df = 0.0
        @constraint model sθ1[N] == sin(θ1f)
        @constraint model cθ1[N] == cos(θ1f)
        @constraint model θ1d[N] == θ1df
        @constraint model sθ2[N] == sin(θ2f)
        @constraint model cθ2[N] == cos(θ2f)
        @constraint model θ2d[N] == θ2df 

        # Objective
        if objectivetype isa MinEffort
            @objective model Min τ ⋅ τ
        elseif objectivetype isa MinTime
            @objective model Min sum(Δt)
        else
            throw(ArgumentError("Objective type not recognized"))
        end

        new(parameters, model, Δt,  sΔθ1,  cΔθ1,  Δθ1, sθ1, cθ1, θ1d, θ1dd, sΔθ2, cΔθ2, Δθ2, sθ2, cθ2, θ2d, θ2dd, τ)
    end
end

end