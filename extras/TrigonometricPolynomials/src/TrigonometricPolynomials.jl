module TrigonometricPolynomials

export
    TrigPoly,
    SinCosDict,
    simplify

using LinearAlgebra
using Reexport
@reexport using DynamicPolynomials

const Var = PolyVar{true}
const Poly{T} = Polynomial{true, T}

struct SinCosVars
    s::Var
    c::Var
end

function SinCosVars(x::Var)
    basename = name(x)
    SinCosVars(Var("sin($basename)"), Var("cos($basename)"))
end

const SinCosDict = Dict{Var, SinCosVars}

struct TrigPoly{T<:Real} <: Real
    poly::Poly{T}
    sincosmap::SinCosDict

    TrigPoly{T}(x::Poly{T}, sincosmap::SinCosDict) where {T<:Real} = new{T}(x, sincosmap)
end

TrigPoly(x::Poly{T}, sincosmap::SinCosDict) where {T<:Real} = TrigPoly{T}(x, sincosmap)
TrigPoly(x, sincosmap::SinCosDict) = TrigPoly(Poly(x), sincosmap)
TrigPoly(x) = TrigPoly(x, SinCosDict())
TrigPoly{T}(x, sincosmap::SinCosDict) where {T<:Real} = TrigPoly(Poly{T}(x), sincosmap)
TrigPoly{T}(x) where {T<:Real} = TrigPoly{T}(x, SinCosDict())
TrigPoly{T}(x::TrigPoly) where {T} = TrigPoly{T}(x.poly, x.sincosmap)

# Utility
function combine_sin_cos_maps(x::SinCosDict, y::SinCosDict)
    if x === y
        x
    elseif isempty(x)
        y
    elseif isempty(y)
        x
    else
        check = (a, b) -> a == b || throw(ArgumentError("Multiple sin/cos variables associated with a variable."))
        merge(check, x, y)
    end
end

function simplify(p::TrigPoly)
    # sin^2 + cos^2 == 1
    poly = p.poly
    for (x, sc) in p.sincosmap
        s = sc.s
        c = sc.c
        if s in variables(poly) && c in variables(poly)
            d, r = divrem(poly, s^2 + c^2)
            poly = d + r
        end
    end
    TrigPoly(poly, p.sincosmap)
end

# Pretty-printing
Base.show(io::IO, p::TrigPoly) = show(io, p.poly)
Base.show(io::IO, mime::MIME"text/plain", p::TrigPoly) = show(io, mime, p.poly)
Base.show(io::IO, mime::MIME"text/latex", p::TrigPoly) = show(io, mime, p.poly)
Base.show(io::IO, mime::MIME"text/print", p::TrigPoly) = show(io, mime, p.poly)

# Math
for (fun, varname) in [(:sin, :s), (:cos, :c)]
    @eval function Base.$fun(p::TrigPoly{T}) where {T<:Real}
        poly = p.poly
        maxdegree(poly) > 1 && throw(ArgumentError("Cannot handle degree > 1."))

        # base case (including constants)
        if nterms(poly) <= 1
            t = term(poly)
            coeff = coefficient(t)
            if isconstant(t)
                return TrigPoly($fun(coeff), p.sincosmap)
            elseif abs(coefficient(t)) == 1
                x = variable(monomial(t))
                sc = get!(() -> SinCosVars(x), p.sincosmap, x)
                R = typeof(sin(coeff))
                if $fun === sin
                    # sin(-x) == -sin(x)
                    return TrigPoly{R}(coeff * sc.$varname, p.sincosmap)
                else
                    # cos(-x) == cos(x)
                    return TrigPoly{R}(abs(coeff) * sc.$varname, p.sincosmap)
                end
            end
        end

        # Recursion, both for handling multiple terms and terms with integer coefficient > 1
        # is handled using the angle-sum formula.
        head = leadingterm(poly)
        tail = removeleadingterm(poly)
        c = coefficient(head)
        isinteger(c) || throw(ArgumentError("Cannot handle non-integer coefficients."))
        if c > 1
            tail += head - monomial(head)
            head = typeof(head)(monomial(head))
        end
        α = TrigPoly(head, p.sincosmap)
        β = TrigPoly(tail, p.sincosmap)
        if $fun === sin
            return sin(α) * cos(β) + cos(α) * sin(β)
        else
            return cos(α) * cos(β) - sin(α) * sin(β)
        end
    end
end

Base.sincos(p::TrigPoly) = (sin(p), cos(p))

for op in [:+, :-, :*]
    @eval Base.$op(x::TrigPoly, y::TrigPoly) = TrigPoly($op(x.poly, y.poly), combine_sin_cos_maps(x.sincosmap, y.sincosmap))
end

for op in [:+, :-, :*, :/, :\]
    @eval Base.$op(x::TrigPoly, y::Real) = TrigPoly($op(x.poly, y), x.sincosmap)
    @eval Base.$op(x::Real, y::TrigPoly) = TrigPoly($op(x, y.poly), y.sincosmap)
end

for op in [:+, :-]
    @eval Base.$op(p::TrigPoly) = TrigPoly($op(p.poly), p.sincosmap)
end

for fun in [:zero, :one]
    @eval Base.$fun(::Type{TrigPoly{T}}) where {T<:Real} = TrigPoly($fun(T))
    @eval Base.$fun(p::TrigPoly) = TrigPoly($fun(p.poly), p.sincosmap)
end

Base.:^(x::TrigPoly, p::Integer) = Base.power_by_squaring(x, p)


# Number-like interface
LinearAlgebra.dot(x::TrigPoly, y::TrigPoly) = x * y
LinearAlgebra.dot(x::TrigPoly, y::Number) = x * y
LinearAlgebra.dot(x::Number, y::TrigPoly) = x * y
LinearAlgebra.symmetric_type(::Type{T}) where {T<:TrigPoly} = T
Base.transpose(p::TrigPoly) = TrigPoly(transpose(p.poly), p.sincosmap)
Base.adjoint(p::TrigPoly) = TrigPoly(adjoint(p.poly), p.sincosmap)
Base.broadcastable(p::TrigPoly) = Ref(p)
Base.conj(p::TrigPoly) = p

# Promotion/conversion
Base.promote_rule(::Type{TrigPoly{T1}}, ::Type{TrigPoly{T2}}) where {T1<:Real, T2<:Real} = TrigPoly{promote_type(T1, T2)}
Base.promote_rule(::Type{TrigPoly{T1}}, ::Type{T2}) where {T1<:Real, T2<:Real} = TrigPoly{promote_type(T1, T2)}

Base.convert(::Type{TrigPoly{T}}, x::TrigPoly{T}) where {T<:Real} = x
Base.convert(::Type{TrigPoly{T}}, x::TrigPoly) where {T<:Real} = TrigPoly{T}(x)
Base.convert(::Type{TrigPoly{T}}, x::Real) where {T} = TrigPoly{T}(x)

end # module
