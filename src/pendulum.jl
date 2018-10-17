@with_kw struct PendulumParameters
    m::Float64 = 1.0
    l::Float64 = 2.0
    g::Float64 = 9.81
    b::Float64 = 0.0#0.1
end

struct PendulumSwingUpProblem
    parameters::PendulumParameters
    model::JuMP.Model
    sΔθ::Vector{JuMP.Variable}
    cΔθ::Vector{JuMP.Variable}
    Δθ::Vector{JuMP.Variable}
    sθ::Vector{JuMP.Variable}
    cθ::Vector{JuMP.Variable}
    θd::Vector{JuMP.GenericAffExpr{Float64,Variable}}
    θdd::Vector{JuMP.Variable}
    τ::Vector{JuMP.Variable}

    function PendulumSwingUpProblem(parameters::PendulumParameters, x0, solver;
            τmax::Number, # maximum torque
            N::Integer, # number of integration steps
            Δt::Number, # time step
            Δθmax::Number = 0.5
        )
        model = Model(solver=solver)

        # Initial state
        θ0, θd0 = x0
        sθprev, cθprev = sincos(θ0)
        θdprev = θd0

        # Variable / expression vectors
        sΔθ = Vector{JuMP.Variable}(undef, N)
        cΔθ = Vector{JuMP.Variable}(undef, N)
        Δθ = Vector{JuMP.Variable}(undef, N)
        sθ = Vector{JuMP.Variable}(undef, N)
        cθ = Vector{JuMP.Variable}(undef, N)
        θd = Vector{JuMP.GenericAffExpr{Float64,Variable}}(undef, N)
        θdd = Vector{JuMP.Variable}(undef, N)
        τ = Vector{JuMP.Variable}(undef, N)

        # Add stages
        for i = 1 : N
            # Kinematics delta
            sΔθ[i], cΔθ[i] = sincosvar(model, "Δθ_$i", θmax=Δθmax, normconstraint=false)
            Δθ[i] = sΔθ[i] # first-order approximation

            # Absolute kinematics (from rotation matrix product)
            sθ[i], cθ[i] = sincosvar(model, "θ_$i", normconstraint=true)
            @NLconstraint model sθ[i] == sθprev * cΔθ[i] + cθprev * sΔθ[i]
            @NLconstraint model cθ[i] == cθprev * cΔθ[i] - sθprev * sΔθ[i]
            sθprev, cθprev = sθ[i], cθ[i]

            # Velocity, acceleration
            θd[i] = @expression model Δθ[i] / Δt
            θdd[i] = @variable model basename="θdd_$i"
            @constraint model θdd[i] * Δt == θd[i] - θdprev
            θdprev = θd[i]

            # Torque
            τ[i] = @variable model basename="τ_$i"
            setlowerbound(τ[i], -τmax)
            setupperbound(τ[i],  τmax)

            # Dynamics
            M = parameters.m * parameters.l^2
            mgl = parameters.m * parameters.g * parameters.l
            b_over_Δt = parameters.b / Δt
            c = @NLexpression model b_over_Δt * Δθ[i] + mgl * sθ[i] # want to use b * θdi, but JuMP won't let me.
            @NLconstraint model M * θdd[i] + c == τ[i]
        end

        # Final state constraint
        θf = π
        θdf = 0.0
        @constraint model sθ[N] == sin(θf)
        @constraint model cθ[N] == cos(θf)
        @constraint model θd[N] == θdf

        # Set objective
        @objective model Min τ ⋅ τ

        new(parameters, model, sΔθ, cΔθ, Δθ, sθ, cθ, θd, θdd, τ)
    end
end
