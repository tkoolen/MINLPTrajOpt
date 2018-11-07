module TrigonometricPolynomialTest

using TrigonometricPolynomials
using Test

@testset "basics" begin
    @polyvar x
    p = TrigPoly(x)
    @test p.poly == polynomial(x)
    sx = sin(p)
    @test haskey(p.sincosmap, x)
    @test name(p.sincosmap[x].s) == "sin(x)"
    @test name(p.sincosmap[x].c) == "cos(x)"
    @test sx.poly == p.sincosmap[x].s
end

@testset "angle sum" begin
    @polyvar α β
    p = TrigPoly(α + β)
    sp = sin(p)
    cp = cos(p)
    sα, cα = p.sincosmap[α].s, p.sincosmap[α].c
    sβ, cβ = p.sincosmap[β].s, p.sincosmap[β].c
    @test sp.poly == sα * cβ + cα * sβ
    @test cp.poly == cα * cβ - sα * sβ
end

@testset "reflection" begin
    @polyvar θ
    p = TrigPoly(-θ)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = p.sincosmap[θ].s, p.sincosmap[θ].c
    @test sp.poly == -sθ
    @test cp.poly == cθ
end

@testset "double angle" begin
    @polyvar θ
    p = TrigPoly(2θ)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = p.sincosmap[θ].s, p.sincosmap[θ].c
    @test sp.poly == 2 * sθ * cθ
    @test cp.poly == cθ^2 - sθ^2
end

@testset "triple angle" begin
    @polyvar θ
    p = TrigPoly(3θ)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = p.sincosmap[θ].s, p.sincosmap[θ].c
    @test sp.poly == 3 * sθ * cθ^2 - sθ^3 # == 3 * sθ - 4 * cθ
    @test cp.poly == cθ^3 - 3 * sθ^2 * cθ
end

@testset "show" begin
    @polyvar θ
    p = sin(TrigPoly(3θ))
    buf = IOBuffer()
    show(buf, p)
    result = String(take!(buf))
    @test occursin("-sin(θ)³ + 3.0sin(θ)cos(θ)²", result)
end

@testset "combine" begin
    @polyvar x y
    p1 = sin(TrigPoly(x))
    p2 = cos(TrigPoly(y + 3))
    p3 = p1 + p2
    @test p3.sincosmap[x] === p1.sincosmap[x]
    @test p3.sincosmap[y] === p2.sincosmap[y]

    p4 = cos(TrigPoly(x))
    @test_throws ArgumentError p1 + p4

    @test (p3 + zero(TrigPoly{Int})).poly == p3.poly
end

end # module
