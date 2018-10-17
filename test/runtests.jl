using MINLPTrajOpt
using Test
using JuMP
using BARON

@testset "pendulum" begin
    N = 30
    θ0 = -0.1
    θd0 = 0.0
    T = 4.5
    Δt = T / N
    parameters = PendulumParameters()
    τmax = parameters.m * parameters.g * parameters.l * 0.9
    solver = BaronSolver(threads=Sys.CPU_THREADS)
    problem = PendulumSwingUpProblem(parameters, (θ0, θd0), solver;
        τmax = τmax, Δt=Δt, N=N)
    solve(problem.model)

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
