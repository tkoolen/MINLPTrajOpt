function sincosvar(model::JuMP.Model, basename::AbstractString; normconstraint=true, θmax=Inf)
    s = @variable model basename="s"*basename
    c = @variable model basename="c"*basename
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
    s, c
end
