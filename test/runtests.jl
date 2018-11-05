module PendulumTest

using MINLPTrajOpt
using MINLPTrajOpt.Pendulum
using Test
using JuMP
using BARON

@testset "pendulum" begin
    solver = BaronSolver(threads=Sys.CPU_THREADS)
    parameters = PendulumParameters()
    τmax = parameters.m * parameters.g * parameters.l * 0.7
    θ0 = 0.1
    θd0 = 0.0

    N = 30
    objectivetype = MinEffort()
    T = 4.0
    Δt = T / N
    Δtmin = Δt
    Δtmax = Δt

    problem = PendulumSwingUpProblem(parameters, (θ0, θd0), solver;
    τmax=τmax, Δtmin=Δtmin, Δtmax=Δtmax, T=T, N=N, objectivetype=objectivetype)
    status = solve(problem.model)

    @test status == :Optimal

    for i = 1 : N
        @test getvalue(problem.sΔθ[i])^2 + getvalue(problem.cΔθ[i])^2 ≈ 1 atol=1e-5
        @test getvalue(problem.sθ[i])^2 + getvalue(problem.cθ[i])^2 ≈ 1 atol=1e-5
    end
    @test diff(getvalue.(problem.θd)) ./ Δt - getvalue(problem.θdd)[2 : end] ≈ zeros(N - 1) atol=1e-6

    Δθ = atan.(getvalue.(problem.sΔθ) ./ getvalue.(problem.cΔθ))
    θ = θ0 .+ [0.; cumsum(Δθ)]
    @test cos(θ[end]) ≈ cos(π) atol=1e-6
    @test sin(θ[end]) ≈ sin(π) atol=1e-6
end

end # module PendulumTest

let
    notebookdir = joinpath(@__DIR__, "..", "notebooks")
    excludedirs = [".ipynb_checkpoints"]
    excludefiles = String[]
    for (root, dir, files) in walkdir(notebookdir)
        basename(root) in excludedirs && continue
        for file in files
            file in excludefiles && continue
            name, ext = splitext(file)
            lowercase(ext) == ".ipynb" || continue
            path = joinpath(root, file)
            @eval module $(gensym()) # Each notebook is run in its own module.
            using Test
            using NBInclude
            @testset "Notebook: $($name)" begin
                # Note: use #NBSKIP in a cell to skip it during tests.
                @nbinclude($path; regex = r"^((?!\#NBSKIP).)*$"s)
            end
            end # module
        end
    end
end
