struct MinEffort <: AbstractObjective end

function setobjective!(problem::TrajOptProblem, ::MinEffort)
    τs = problem.τs
    N = length(τs)
    @objective problem.model Min sum(τs[i] ⋅ τs[i] for i = 1 : N) / N
end

