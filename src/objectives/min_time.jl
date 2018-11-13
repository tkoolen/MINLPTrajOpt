struct MinTime <: AbstractObjective end

function setobjective!(problem::TrajOptProblem, ::MinTime)
    @objective problem.model Min sum(problem.Δts)
end

