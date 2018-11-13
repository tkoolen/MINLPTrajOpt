function sincosconstraints(model::Model, s::Variable, c::Variable; normconstraint=true, θmax=Inf)
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

