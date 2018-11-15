function sincosconstraints(model::Model, s::Variable, c::Variable; normconstraint=true, θmax=Inf)
    # Add bounds to help the solver out:
    setlowerbound.((s, c), -1.0)
    setupperbound.((s, c),  1.0)

    if normconstraint
        @NLconstraint model s^2 + c^2 == 1
    end

    if θmax < π / 2
        @constraint model -sin(θmax) <= s <= sin(θmax)
        if θmax < π
            @constraint model cos(θmax) <= c
        end
    end
    nothing
end

function quatconstraints(model::JuMP.Model, q::Quat{Variable}; normconstraint=true, θmax=π, w_nonnegative=false)
    w = q.w
    x = q.x
    y = q.y
    z = q.z

    # Add bounds to help the solver out:
    setlowerbound.((w, x, y, z), -1.0)
    setupperbound.((w, x, y, z),  1.0)

    if w_nonnegative
        # Avoid unnecessary branches due to Quat double cover.
        setlowerbound(w, 0.0)
    end

    if normconstraint
        @NLconstraint model w^2 + x^2 + y^2 + z^2 == 1
    end

    if θmax < π / 2
        # TODO: norm constraint for the vector part?
        setlowerbound.((x, y, z), -sin(θmax / 2))
        setupperbound.((x, y, z),  sin(θmax / 2))
        if θmax < π
            setlowerbound(w, cos(θmax / 2))
        end
    end
    nothing
end

function quatmulconstraints(model::JuMP.Model, q1::Quat, q2::Quat, q12::Quat)
    @NLconstraint model q12.w == q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
    @NLconstraint model q12.x == q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y
    @NLconstraint model q12.y == q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x
    @NLconstraint model q12.z == q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    nothing
end

