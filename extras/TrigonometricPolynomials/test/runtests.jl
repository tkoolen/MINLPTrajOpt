module TrigonometricPolynomialTest

using TrigonometricPolynomials
using Test

@testset "basics" begin
    @polyvar x
    sincosmap = SinCosDict()
    p = TrigPoly(x, sincosmap)
    @test p.poly == polynomial(x)
    sx = sin(p)
    @test haskey(sincosmap, x)
    @test name(sincosmap[x].s) == "sin(x)"
    @test name(sincosmap[x].c) == "cos(x)"
    @test sx.poly == sincosmap[x].s
end

@testset "angle sum" begin
    @polyvar α β
    sincosmap = SinCosDict()
    p = TrigPoly(α + β, sincosmap)
    sp = sin(p)
    cp = cos(p)
    sα, cα = sincosmap[α].s, sincosmap[α].c
    sβ, cβ = sincosmap[β].s, sincosmap[β].c
    @test sp.poly == sα * cβ + cα * sβ
    @test cp.poly == cα * cβ - sα * sβ
end

@testset "reflection" begin
    @polyvar θ
    sincosmap = SinCosDict()
    p = TrigPoly(-θ, sincosmap)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = sincosmap[θ].s, sincosmap[θ].c
    @test sp.poly == -sθ
    @test cp.poly == cθ
end

@testset "double angle" begin
    @polyvar θ
    sincosmap = SinCosDict()
    p = TrigPoly(2θ, sincosmap)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = sincosmap[θ].s, sincosmap[θ].c
    @test sp.poly == 2 * sθ * cθ
    @test cp.poly == cθ^2 - sθ^2
end

@testset "triple angle" begin
    @polyvar θ
    sincosmap = SinCosDict()
    p = TrigPoly(3θ, sincosmap)
    sp = sin(p)
    cp = cos(p)
    sθ, cθ = sincosmap[θ].s, sincosmap[θ].c
    @test sp.poly == 3 * sθ * cθ^2 - sθ^3 # == 3 * sθ - 4 * cθ
    @test cp.poly == cθ^3 - 3 * sθ^2 * cθ
end

@testset "show" begin
    @polyvar θ
    sincosmap = SinCosDict()
    p = sin(TrigPoly(3θ, sincosmap))
    buf = IOBuffer()
    show(buf, p)
    @test occursin(String(take!(buf)), "TrigPoly{Float64}: -sin(θ)³ + 3.0sin(θ)cos(θ)²")
end

end # module
