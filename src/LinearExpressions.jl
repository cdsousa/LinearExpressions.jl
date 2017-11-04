module LinearExpressions


export Symbolic, LinExpr, AbstractVariable, BasicVariable, TypedVariable,
    RealVariable, getcoeff, setcoeff!, differentiate


import Base: hash, show, print, showcompact, convert, promote_rule, zero, one,
             conj, (==), (!=), (+), (-), (*), (.*), transpose



# Borrowed from Calculus.jl
abstract type Symbolic end
abstract type AbstractVariable <: Symbolic end
struct BasicVariable <: AbstractVariable
  sym::Symbol
end
show(io::IO, x::BasicVariable) = print(io, x.sym)
(==)(x::BasicVariable, y::BasicVariable) = x.sym == y.sym


# Symbolic variables with assosiated type (experimental)
struct TypedVariable{T} <: AbstractVariable
  sym::Symbol
end
show(io::IO, x::TypedVariable) = print(io, x.sym)
(==)(x::TypedVariable{T}, y::TypedVariable{T}) where {T} = x.sym == y.sym
const RealVariable = TypedVariable{Real}


###########################################################


const ScalarCoeff = Real
const ArrayCoeff{T <: Real} = AbstractArray{T}
const Coeff = Union{ScalarCoeff,ArrayCoeff}



mutable struct LinExpr{Tc<:Coeff,Tv<:AbstractVariable}
  constt::Tc
  coeffs::Dict{Tv,Tc}

  LinExpr{Tc,Tv}(c::Tc, vc::Dict{Tv,Tc}) where {Tc,Tv} = new(c, filter((v, c) -> c != zero(c), vc))
  LinExpr{Tc,Tv}(c::Tc) where {Tc,Tv} = new(c, Dict{Tv,Tc}())
  LinExpr{Tc,Tv}() where {Tc,Tv} = new(zero(Tc), Dict{Tv,Tc}())
end

LinExpr{Tc,Tv}(vc::Dict{Tv,Tc}) where {Tc,Tv} = LinExpr{Tc,Tv}(zero(Tc), vc)
LinExpr(c::Tc, vc::Dict{Tv,Tc}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(c, vc)
LinExpr(vc::Dict{Tv,Tc}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(zero(Tc), vc)



function linexpr_show(io::IO, e::LinExpr{Tc,Tv}, compact::Bool) where {Tc <: ScalarCoeff, Tv <: AbstractVariable}
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
show(io::IO, e::LinExpr{Tc,Tv}) where {Tc <: ScalarCoeff, Tv <: AbstractVariable} = linexpr_show(io, e, false)
showcompact(io::IO, e::LinExpr{Tc,Tv}) where {Tc <: ScalarCoeff, Tv <: AbstractVariable} = linexpr_show(io, e, true)


function show(io::IO, e::LinExpr{Tc,Tv}) where {Tc <: ArrayCoeff, Tv <: AbstractVariable}
  start = true
  if e.constt != zero(e.constt) || isempty(e.coeffs)
    print(io, e.constt)
    start = false
  end
  for (s, c) in e.coeffs
    if start
      start = false
    else
      print(io, "\n+\n\n")
    end
    print(io, c)
    print(io, s)
    print(io, "\n")
  end
end


function ==(a::LinExpr, b::LinExpr)
  a.constt == b.constt && a.coeffs == b.coeffs
end


function getcoeff(e::LinExpr{Tc,Tv}, v::Tv) where {Tc <: Coeff, Tv <: AbstractVariable}
    get(e.coeffs, v, zero(e.constt))
end
function setcoeff!(e::LinExpr{Tc,Tv}, v::Tv, x::Tc) where {Tc <: Coeff, Tv <: AbstractVariable}
    e.coeffs[v] = x
end


one(e::LinExpr{Tc,Tv}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(one(e.constt))
zero(e::LinExpr{Tc,Tv}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(zero(e.constt))

one(::Type{LinExpr{Tc,Tv}}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(one(Tc))
zero(::Type{LinExpr{Tc,Tv}}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(zero(Tc))

conj(x::TypedVariable{T}) where {T <: Real} = x
conj(e::LinExpr{Tc,TypedVariable{T}}) where {Tc <: Real, T <: Real} = e


function convert(::Type{Tc}, e::LinExpr{Tc,Tv}) where {Tc <: Coeff, Tv <: AbstractVariable}
    if isempty(e.coeffs)
        e.constt
    else
        throw(InexactError())
    end
end

function convert(::Type{Tv}, e::LinExpr{Tc,Tv}) where {Tc <: ScalarCoeff, Tv <: AbstractVariable}
    if e.constt != zero(Tc) || length(e.coeffs) != 1
        throw(InexactError())
    else
        s = collect(keys(e.coeffs))[1]
        if e.coeffs[s] != one(Tc)
            throw(InexactError())
        else
            s
        end
    end
end



convert(::Type{LinExpr{Tc,Tv}}, e::LinExpr{Tc,Tv}) where {Tc <: Coeff, Tv <: AbstractVariable} = e
function convert(::Type{LinExpr{Tc2,Tv}}, e::LinExpr{Tc1,Tv}) where {Tc1 <: Coeff, Tc2 <: Coeff, Tv <: AbstractVariable}
    LinExpr{Tc2,Tv}(convert(Tc2, e.constt), Dict(v => convert(Tc2, c) for (v, c) in e.coeffs))
end


convert(::Type{LinExpr{Tc1,Tv}}, x::Tc2) where {Tc1 <: Coeff, Tv <: AbstractVariable, Tc2 <: Coeff} = LinExpr{Tc1,Tv}(convert(Tc1, x), Dict{Tv,Tc1}())
convert(::Type{LinExpr{Tc,Tv}}, x::Tv) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(zero(Tc), Dict(x => one(Tc)))


promote_rule(::Type{LinExpr{Tc1,Tv}}, ::Type{LinExpr{Tc2,Tv}}) where {Tc1 <: Coeff, Tc2 <: Coeff, Tv <: AbstractVariable} = LinExpr{promote_type(Tc1, Tc2),Tv}
promote_rule(::Type{LinExpr{Tc1,Tv}}, ::Type{Tc2}) where {Tc1 <: Coeff, Tv <: AbstractVariable, Tc2 <: Coeff} = LinExpr{promote_type(Tc1, Tc2),Tv}
promote_rule(::Type{Tc}, ::Type{Tv}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}
promote_rule(::Type{LinExpr{Tc,Tv}}, ::Type{Tv}) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}


function promote_rule(::Type{LinExpr{Array{T1,N},Tv}}, ::Type{LinExpr{Array{T2,N},Tv}}) where {T1 <: ScalarCoeff, T2 <: ScalarCoeff, N, Tv <: AbstractVariable}
  LinExpr{Array{promote_type(T1, T2),N},Tv}
end
function promote_rule(::Type{LinExpr{Array{T1,N},Tv}}, ::Type{Array{T2,N}}) where {T1 <: ScalarCoeff, T2 <: ScalarCoeff, N, Tv <: AbstractVariable}
  LinExpr{Array{promote_type(T1, T2),N},Tv}
end


function *(c::Tc, e::LinExpr{Tc,Tv}) where {Tc <: Coeff, Tv <: AbstractVariable}
  constt = e.constt * c
  coeffs = similar(e.coeffs)
  for (s, coeff) in e.coeffs
    coeffs[s] = c * coeff
  end
  LinExpr{Tc,Tv}(constt, coeffs)
end
*(e::LinExpr{Tc,Tv}, c::Tc) where {Tc <: Coeff, Tv <: AbstractVariable} = c * e

function *(c::Tc2, e::LinExpr{Tc1,Tv}) where {Tc1 <: Coeff, Tv <: AbstractVariable, Tc2 <: Coeff}
    Tc_common = promote_type(Tc1, Tc2)
    convert(Tc_common, c) * convert(LinExpr{Tc_common,Tv}, e)
end
*(e::LinExpr{Tc1,Tv}, c::Tc2) where {Tc1 <: Coeff, Tv <: AbstractVariable, Tc2 <: Coeff} = c * e

*(c::Tc, s::Tv) where {Tc <: Coeff, Tv <: AbstractVariable} = LinExpr{Tc,Tv}(zero(c), Dict(s => c))
*(s::Tv, c::Tc) where {Tc <: Coeff, Tv <: AbstractVariable} = c * s




const VarOrLinExpr = Union{AbstractVariable,LinExpr}



# .*(a::ScalarCoeff, e::VarOrLinExpr) = a * e
# .*(e::VarOrLinExpr, a::ScalarCoeff) = e * a
#
# .*(e::VarOrLinExpr, a::ArrayCoeff) = map(x -> e * x, a)
# .*(a::ArrayCoeff, e::VarOrLinExpr) = e .* a



function +(a::T, b::T) where {T <: LinExpr}
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

-(a::T, b::T) where {T <: LinExpr} = a + (-b)

-(x::AbstractVariable) = -1 * x
+(x::AbstractVariable) = 1 * x


+(x::VarOrLinExpr, y::VarOrLinExpr) = +(promote(+x, y)...)
-(x::VarOrLinExpr, y::VarOrLinExpr) = x + (-y)

+(x::VarOrLinExpr, y::Coeff) = +(promote(x, y)...)
+(x::Coeff, y::VarOrLinExpr) = y + x

-(x::VarOrLinExpr, y::Coeff) = -(promote(x, y)...)
-(x::Coeff, y::VarOrLinExpr) = x + (-y)


transpose(x::RealVariable) = x
transpose(x::LinExpr{<:Real, <:RealVariable}) = x


==(x::VarOrLinExpr, y::Coeff) = ==(promote(x, y)...)
!=(x::VarOrLinExpr, y::Coeff) = !(==(promote(x, y)...))


==(::T1, ::T2) where {T1 <: AbstractVariable, T2 <: AbstractVariable} = false  # safe for different variable types
==(x::VarOrLinExpr, y::VarOrLinExpr) = ==(promote(x, y)...)
!=(x::VarOrLinExpr, y::VarOrLinExpr) = !(==(promote(x, y)...))


differentiate(s::AbstractVariable, d::AbstractVariable) = s == d ? one(1) : zero(0)
differentiate(s::TypedVariable{T}, d::TypedVariable{T}) where {T} = s == d ? one(T) : zero(T)
differentiate(c::Coeff, d::AbstractVariable) = zero(c)
differentiate(e::LinExpr{Tc,Tv}, d::Tv) where {Tc, Tv} = get(e.coeffs, d, zero(e.constt))


end # module
