abstract type AbstractObjective end
function setobjective! end

struct TrajOptProblem
    mechanism::Mechanism{Float64}
    model::JuMP.Model
    Δts::Vector{JuMP.Variable}
    Δqs::Vector{JointSegmentedVector{JuMP.Variable}}
    qs::Vector{JointSegmentedVector{JuMP.Variable}}
    vs::Vector{JointSegmentedVector{JuMP.Variable}}
    v̇s::Vector{JointSegmentedVector{JuMP.Variable}}
    τs::Vector{JointSegmentedVector{JuMP.Variable}}

    function TrajOptProblem(mechanism::Mechanism, x0::MechanismState, solver;
            xf::Union{MechanismState, Nothing}=nothing, # final state
            N::Integer, # number of integration steps
            Δtmin::Number, # min time step
            Δtmax::Number, # max time step
            objective::AbstractObjective,
            T::Union{Number, Nothing} = objectivetype isa MinTime ? nothing : N * Δtmax, # final time
            Δθmax::Number = 0.5, # max angle change during one integration time step (to control badness of small-angle approximation)
        )
        model = Model(solver=solver)
        nq = num_positions(mechanism)
        nv = num_velocities(mechanism)

        # Time step settings
        @assert Δtmax >= Δtmin
        fixedstep = Δtmin == Δtmax
        if fixedstep
            T == N * Δtmax || throw(ArgumentError("Total time does not match Δt."))
        end

        # Initial state
        qprev = configuration(x0)
        vprev = velocity(x0)

        # Variable vectors
        Δts = Vector{JuMP.Variable}(undef, N)
        Δqs = [similar(qprev, JuMP.Variable) for i = 1 : N]
        qs = [similar(qprev, JuMP.Variable) for i = 1 : N]
        vs = [similar(vprev, JuMP.Variable) for i = 1 : N]
        v̇s = [similar(vprev, JuMP.Variable) for i = 1 : N]
        τs = [similar(vprev, JuMP.Variable) for i = 1 : N]

        # Symbolic dynamics
        qsym, vsym, v̇sym, τsym, dynsym = symbolic_dynamics(mechanism)

        for i = 1 : N
            # Create variables
            Δts[i] = @variable model basename="Δt_{$i}"
            Δq = Δqs[i] .= [@variable model basename="Δq_{$j,$i}" for j = 1 : nq]
            q = qs[i] .= [@variable model basename="q_{$j,$i}" for j = 1 : nq]
            v = vs[i] .= [@variable model basename="v_{$j,$i}" for j = 1 : nv]
            v̇ = v̇s[i] .= [@variable model basename="v̇_{$j,$i}" for j = 1 : nv]
            τ = τs[i] .= [@variable model basename="τ_{$j,$i}" for j = 1 : nv]

            # Time steps
            if fixedstep
                JuMP.fix(Δts[i], Δtmin)
                Δt = Δtmin
            else
                JuMP.setlowerbound(Δts[i], Δtmin)
                JuMP.setupperbound(Δts[i], Δtmax)
                Δt = Δts[i]
            end

            for joint in tree_joints(mechanism)
                if joint_type(joint) isa SinCosRevolute
                    # Kinematics delta
                    sΔθ, cΔθ = Δq[joint]
                    sincosconstraints(model, sΔθ, cΔθ; normconstraint=false, θmax=Δθmax)

                    # Absolute kinematics (from rotation matrix product)
                    sθprev, cθprev = qprev[joint]
                    sθ, cθ = q[joint]
                    sincosconstraints(model, sθ, cθ; normconstraint=true)
                    @NLconstraint model sθ == sθprev * cΔθ + cθprev * sΔθ
                    @NLconstraint model cθ == cθprev * cΔθ - sθprev * sΔθ

                    # Velocity
                    Δθ = sΔθ # first-order approximation
                    θd = v[joint][1]
                    @constraint model Δt * θd == Δθ
                elseif joint_type(joint) isa QuaternionSpherical
                    # Kinematics delta
                    Δquat = Quat(Δq[joint]..., false)
                    quatconstraints(model, Δquat; normconstraint=false, θmax=Δθmax, w_nonnegative=true)

                    # Absolute kinematics
                    quatprev = Quat(qprev[joint]..., false)
                    quat = Quat(q[joint]..., false)
                    quatconstraints(model, quat; normconstraint=true)
                    quatmulconstraints(model, quatprev, Δquat, quat)

                    # Velocity
                    ω = v[joint]
                    @constraint model Δt .* ω .== [Δquat.x, Δquat.y, Δquat.z] # first-order approximation
                else
                    throw(ArgumentError("Unhandled joint type"))
                end

                # Effort bounds
                foreach(τ[joint], effort_bounds(joint)) do τ, bounds
                    if lower(bounds) == upper(bounds)
                        JuMP.fix(τ, lower(bounds))
                    else
                        setlowerbound(τ, lower(bounds))
                        setupperbound(τ, upper(bounds))
                    end
                end
            end

            # Acceleration
            @constraint model Δt .* v̇ .== v .- vprev

            # Dynamics
            varmap = Dict(Iterators.flatten((
                    zip(variable.(qsym), q),
                    zip(variable.(vsym), v),
                    zip(variable.(v̇sym), v̇),
                    zip(variable.(τsym), τ)
            )))
            JuMP.addNLconstraint.(Ref(model), dynamics_constraints(dynsym, varmap))

            qprev = q
            vprev = v
        end

        # Total time constraint
        if !fixedstep && T !== nothing
            @constraint model sum(Δts) == T
        end

        # Final state constraint
        if xf !== nothing
            @constraint model qs[end] .== configuration(xf)
            @constraint model vs[end] .== velocity(xf)
        end

        problem = new(mechanism, model, Δts, Δqs, qs, vs, v̇s, τs)

        # Objective
        setobjective!(problem, objective)

        problem
    end
end

function symbolic_dynamics(mechanism::Mechanism{X}, simplify_iters=10) where X # TODO: simplify_iters
    nq = num_positions(mechanism)
    nv = num_velocities(mechanism)
    sincosmap = SinCosDict()
    T = TrigPoly{X}
    state = MechanismState{T}(mechanism)
    q = configuration(state)
    v = velocity(state)
    v̇ = similar(v)
    τ = similar(v)
    for (i, joint) in enumerate(joints(mechanism))
        if joint_type(joint) isa SinCosRevolute
            θjoint = T(Var("θ$i"), sincosmap)
            set_configuration!(state, joint, θjoint)
        else
            qjoint = [T(Var("q$i,$j"), sincosmap) for j = 1 : num_positions(joint)]
            set_configuration!(state, joint, qjoint)
        end
        copyto!(v[joint], [T(Var("v$i,$j"), sincosmap) for j = 1 : num_velocities(joint)])
        copyto!(v̇[joint], [T(Var("v̇$i,$j"), sincosmap) for j = 1 : num_velocities(joint)])
        copyto!(τ[joint], [T(Var("τ$i,$j"), sincosmap) for j = 1 : num_velocities(joint)])
    end
    setdirty!(state)
    M = mass_matrix(state)
    c = dynamics_bias(state)
    eom = M * v̇ + c - τ
    for i = 1 : simplify_iters
        for joint in joints(mechanism)
            if joint_type(joint) isa QuaternionSpherical
                w, x, y, z = map(x -> x.poly, q[joint])
                normsquared = w^2 + x^2 + y^2 + z^2 # == 1
                map!(eom, eom) do trigpoly
                    poly = trigpoly.poly
                    d, r = divrem(poly, normsquared)
                    TrigPoly(d + r, trigpoly.sincosmap)
                end
            end
        end
    end
    q, v, v̇, τ, eom
end

function dynamics_constraints(sym_dynamics::AbstractVector{<:TrigPoly}, varmap::Dict{Var})
    nv = length(sym_dynamics)
    exprs = Vector{Expr}(undef, nv)
    for i in 1 : nv
        sum = Expr(:call, :+)
        for t in terms(sym_dynamics[i].poly)
            c = coefficient(t)
            m = monomial(t)
            prod = Expr(:call, :*, coefficient(t))
            for (x, power) in powers(m)
                if power == 1
                    push!(prod.args, :($(varmap[x])))
                elseif power != 0
                    push!(prod.args, :($(varmap[x])^$(power)))
                end
            end
            push!(sum.args, prod)
        end
        exprs[i] = :($sum == 0)
    end
    exprs
end

