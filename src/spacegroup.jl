import Base.*

# Group action is x ↦ a*x+b
struct SpaceGroupElement{N,T<:Integer} <: GroupElement
    a::SMatrix{N,N,T}
    b::SVector{N, Rational{T}}
    function SpaceGroupElement{N,T}(a::SMatrix{N,N,T}, b::SVector{N,Rational{T}}) where {N, T<:Integer}
        new(a, b-floor.(b))
    end
end

function SpaceGroupElement(a::SMatrix{N,N,T}, b::SVector{N,Rational{T}}) where {N,T<:Integer}
    SpaceGroupElement{N,T}(a, b)
end

SpaceGroupElement{N,T}() where {N,T<:Integer} =
    SpaceGroupElement(one(SMatrix{N,N,T}), zeros(SVector{N, Rational{T}}))

function *(e1::SpaceGroupElement{N,T}, e2::SpaceGroupElement{N,T}) where {N,T<:Integer}
    SpaceGroupElement(e1.a*e2.a, e1.a*e2.b+e1.b)
end

const SpaceGroup{N,T} = FiniteGroup{SpaceGroupElement{N,T}}