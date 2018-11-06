module TrigonometricPolynomials

export
    TrigPoly,
    SinCosDict

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

struct TrigPoly{T<:Real}
    poly::Poly{T}
    sincosmap::SinCosDict

    TrigPoly{T}(x::Poly{T}, sincosmap::SinCosDict) where {T<:Real} = new{T}(x, sincosmap)
end

TrigPoly(x::Poly{T}, sincosmap::SinCosDict) where {T<:Real} = TrigPoly{T}(x, sincosmap)
TrigPoly(x, sincosmap::SinCosDict) = TrigPoly(Poly(x), sincosmap)
TrigPoly{T}(x, sincosmap::SinCosDict) where {T<:Real} = TrigPoly(Poly{T}(x), sincosmap)
TrigPoly{T}(x::TrigPoly) where {T} = TrigPoly{T}(x.poly, x.sincosmap)

function print_trigpoly(io::IO, mime::MIME, p::TrigPoly)
    print(io, typeof(p))
    print(io, ": ")
    show(io, mime, p.poly)
    io
end

Base.show(io::IO, p::TrigPoly) = show(io, MIME"text/plain"(), p)
Base.show(io::IO, mime::MIME"text/plain", p::TrigPoly) = print_trigpoly(io, mime, p)
Base.show(io::IO, mime::MIME"text/latex", p::TrigPoly) = print_trigpoly(io, mime, p)
Base.show(io::IO, mime::MIME"text/print", p::TrigPoly) = print_trigpoly(io, mime, p)

# Math
for (fun, varname) in [(:sin, :s), (:cos, :c)]
    @eval function Base.$fun(p::TrigPoly{T}) where {T}
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

for op in [:+, :-, :*]
    @eval function Base.$op(x::TrigPoly, y::TrigPoly)
        x.sincosmap === y.sincosmap || throw(ArgumentError("sincosmap mismatch"))
        TrigPoly($op(x.poly, y.poly), x.sincosmap)
    end
end

for op in [:+, :-, :*, :/, :\]
    @eval Base.$op(x::TrigPoly, y::Real) = TrigPoly($op(x.poly, y), x.sincosmap)
    @eval Base.$op(x::Real, y::TrigPoly) = TrigPoly($op(x, y.poly), y.sincosmap)
end

for op in [:+, :-]
    @eval Base.$op(p::TrigPoly) = TrigPoly($op(p.poly), p.sincosmap)
end

for fun in [:zero, :one]
    @eval Base.$fun(p::TrigPoly) = TrigPoly($fun(p.poly), p.sincosmap)
end

Base.:^(x::TrigPoly, p::Integer) = Base.power_by_squaring(x, p)

# Promotion/conversion
Base.promote_rule(::Type{TrigPoly{T1}}, ::Type{TrigPoly{T2}}) where {T1<:Real, T2<:Real} = TrigPoly{promote_type(T1, T2)}
Base.promote_rule(::Type{TrigPoly{T1}}, ::Type{T2}) where {T1<:Real, T2<:Real} = TrigPoly{promote_type(T1, T2)}

Base.convert(::Type{TrigPoly{T}}, x::TrigPoly{T}) where {T} = x
Base.convert(::Type{TrigPoly{T}}, x::TrigPoly) where {T} = TrigPoly{T}(x)

end # module
