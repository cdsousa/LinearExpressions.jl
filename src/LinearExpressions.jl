module LinearExpressions

export Symbolic, LinExpr, SymReal

import Base: hash, show, print, showcompact, convert, promote_rule, zero, one


const SYMBOLIC_HASH_SEED = uint(112358)

immutable Symbolic{T}
  sym::Symbol
  _hash::Uint
  function Symbolic(sym::Symbol)
    new(sym, bitmix(bitmix(SYMBOLIC_HASH_SEED, hash(T)), hash(sym)))
  end
end

hash(s::Symbolic) = s._hash

show(io::IO, x::Symbolic) = print(io, x.sym)
showcompact(io::IO, x::Symbolic) = show(io, x)


typealias SymReal Symbolic{Real}


typealias Coeff{T<:Real} Union(T, AbstractArray{T})



type LinExpr{Tc<:Coeff, Tv<:Symbolic}
  constt::Tc
  coeffs::Dict{Tv, Tc}

  LinExpr(c::Tc, vc::Dict{Tv, Tc}) = new(c, filter((v, c) -> c != zero(c), vc))
  LinExpr(c::Tc) = new(c, Dict{Tv, Tc}())
  LinExpr(vc::Dict{Tv, Tc}) = LinExpr{Tc, Tv}(zero(Tc), vc)
  LinExpr() = new(zero(Tc), Dict{Tv, Tc}())
end

LinExpr{Tc<:Coeff, Tv<:Symbolic}(c::Tc, vc::Dict{Tv, Tc}) = LinExpr{Tc, Tv}(c, vc)



function linexpr_show{Tc<:Real, Tv<:Symbolic}(io::IO, e::LinExpr{Tc, Tv}, compact::Bool)
  start = true
  if e.constt != zero(e.constt) || isempty(e.coeffs)
    print(io, e.constt)
    start = false
  end
  for (s, c) in e.coeffs
    if start
      start = false
    else
      if c < zero(c)
        c = -c
        print(io, compact ? "-" : " - ")
      else
        print(io, compact ? "+" : " + ")
      end
    end
    if c == -one(c)
      print(io, "-")
    elseif c != one(c)
      print(io, c)
    end
    print(io, s)
  end
end
show{Tc<:Real, Tv<:Symbolic}(io::IO, e::LinExpr{Tc, Tv}) = linexpr_show(io, e, false)
showcompact{Tc<:Real, Tv<:Symbolic}(io::IO, e::LinExpr{Tc, Tv}) = linexpr_show(io, e, true)


function show{Tc<:AbstractArray, Tv<:Symbolic}(io::IO, e::LinExpr{Tc, Tv})
  start = true
  if e.constt != zero(e.constt) || isempty(e.coeffs)
    print(io, e.constt)
    start = false
  end
  for (s, c) in e.coeffs
    if start
      start = false
    else
      print(io,"\n+\n\n")
    end
    print(io, c)
    print(io, s)
    print(io, "\n")
  end
end


function ==(a::LinExpr, b::LinExpr)
  a.constt == b.constt && a.coeffs == b.coeffs
end


one{Tc<:Coeff, Tv<:Symbolic}(e::LinExpr{Tc, Tv}) = one(e.constt)
zero{Tc<:Coeff, Tv<:Symbolic}(e::LinExpr{Tc, Tv}) = zero(e.constt)

one{Tc<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc, Tv}}) = one(Tc)
zero{Tc<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc, Tv}}) = zero(Tc)


convert{Tc<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc, Tv}}, e::LinExpr{Tc, Tv}) = e
function convert{Tc1<:Coeff, Tc2<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc2, Tv}}, e::LinExpr{Tc1, Tv})
    LinExpr{Tc2, Tv}(convert(Tc2, e.constt), [v=>convert(Tc2, c) for (v, c) in e.coeffs])
end


convert{Tc1<:Coeff, Tv<:Symbolic, Tc2<:Coeff}(::Type{LinExpr{Tc1, Tv}}, x::Tc2) = LinExpr{Tc1, Tv}(convert(Tc1, x), Dict{Tv, Tc1}())
convert{Tc<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc, Tv}}, x::Tv) = LinExpr{Tc, Tv}(zero(Tc), [x=>one(Tc)])


promote_rule{Tc1<:Coeff, Tc2<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc1, Tv}}, ::Type{LinExpr{Tc2, Tv}}) = LinExpr{promote_type(Tc1, Tc2), Tv}
promote_rule{Tc1<:Coeff, Tv<:Symbolic, Tc2<:Coeff}(::Type{LinExpr{Tc1, Tv}}, ::Type{Tc2}) = LinExpr{promote_type(Tc1, Tc2), Tv}
promote_rule{Tc<:Coeff, Tv<:Symbolic}(::Type{Tc}, ::Type{Tv}) = LinExpr{Tc, Tv}
promote_rule{Tc<:Coeff, Tv<:Symbolic}(::Type{LinExpr{Tc, Tv}}, ::Type{Tv}) = LinExpr{Tc, Tv}


function promote_rule{T1<:Real, T2<:Real, N, Tv<:Symbolic}(::Type{LinExpr{Array{T1, N}, Tv}}, ::Type{LinExpr{Array{T2, N}, Tv}})
  LinExpr{Array{promote_type(T1, T2), N}, Tv}
end
function promote_rule{T1<:Real, T2<:Real, N, Tv<:Symbolic}(::Type{LinExpr{Array{T1, N}, Tv}}, ::Type{Array{T2, N}})
  LinExpr{Array{promote_type(T1, T2), N}, Tv}
end


function *{Tc<:Coeff, Tv<:Symbolic}(c::Tc, e::LinExpr{Tc, Tv})
  constt = e.constt * c
  coeffs = similar(e.coeffs)
  for (s, coeff) in e.coeffs
    coeffs[s] = c * coeff
  end
  LinExpr{Tc, Tv}(constt, coeffs)
end
*{Tc<:Coeff, Tv<:Symbolic}(e::LinExpr{Tc, Tv}, c::Tc) = c * e

function *{Tc1<:Coeff, Tv<:Symbolic, Tc2<:Coeff}(c::Tc2, e::LinExpr{Tc1, Tv})
    Tc_common = promote_type(Tc1, Tc2)
    convert(Tc_common, c) * convert(LinExpr{Tc_common, Tv}, e)
end
*{Tc1<:Coeff, Tv<:Symbolic, Tc2<:Coeff}(e::LinExpr{Tc1, Tv}, c::Tc2) = c * e

*{Tc<:Coeff, Tv<:Symbolic}(c::Tc, s::Tv) = LinExpr{Tc, Tv}(zero(c), [s=>c])
*{Tc<:Coeff, Tv<:Symbolic}(s::Tv, c::Tc) = c * s


.*(e::Union(Symbolic, LinExpr), a::AbstractArray) = map(x -> e*x, a)
.*(a::AbstractArray, e::Union(Symbolic, LinExpr)) = e .* a



function +{T<:LinExpr}(a::T, b::T)
  constt = a.constt + b.constt
  coeffs = similar(a.coeffs)
  for (v, c) in a.coeffs
    coeffs[v] = c
  end
  for (v, c) in b.coeffs
    coeffs[v] = c + get(coeffs, v, zero(constt))
  end
  LinExpr(constt, coeffs)
end

+(a::LinExpr) = a

function -(a::LinExpr)
  b = LinExpr(-a.constt, similar(a.coeffs))
  for (v, c) in a.coeffs
    b.coeffs[v] = -c
  end
  b
end

-{T<:LinExpr}(a::T, b::T) = a + (-b)

-(x::Symbolic) = -1 * x
+(x::Symbolic) = 1 * x


+(x::Union(Symbolic, LinExpr), y::Union(Symbolic, LinExpr)) = +(promote(+x,y)...)
-(x::Union(Symbolic, LinExpr), y::Union(Symbolic, LinExpr)) = x + (-y)

+(x::Union(Symbolic, LinExpr), y::Coeff) = +(promote(x,y)...)
+(x::Coeff, y::Union(Symbolic, LinExpr)) = y + x

-(x::Union(Symbolic, LinExpr), y::Coeff) = -(promote(x,y)...)
-(x::Coeff, y::Union(Symbolic, LinExpr)) = x + (-y)


end # module


