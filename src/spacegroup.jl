import Base.*, Base.-


"""
    SpaceGroupElement{N,T<:Integer}

An element of a space group in N dimensions.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

The group action is given by the formula x ↦ a*x+b.

# Fields
- `a::StaticArrays.SMatrix{N,N,T}`: The linear transformation matrix.
- `b::StaticArrays.SVector{N, Rational{T}}`: The translation vector.

The group action is given by the formula x ↦ a*x+b.

# Constructors
- `SpaceGroupElement{N,T}()`: The identity element of the space group.
- `SpaceGroupElement{N,T}(t::SVector{N,T})`: A pure translation element of the space group.
- `SpaceGroupElement{N,T}(m::SMatrix{N,N,T})`: A pure linear transformation element of the space group.

# Example
```julia
julia> using StaticArrays

julia> e1 = SpaceGroupElement{2,Int}(SMatrix{2,2,Int}([1 0; 0 1]))
SpaceGroupElement{2,Int64}(
  [1 0; 0 1],
  [0, 0]
)

julia> e2 = SpaceGroupElement{2,Int}(SVector{2,Int}([1, 1]))
SpaceGroupElement{2,Int64}(
  [1 0; 0 1],
  [1, 1]
)

julia> e1*e2
SpaceGroupElement{2,Int64}(
  [1 0; 0 1],
  [1, 1]
)

julia> e1∘e2
SpaceGroupElement{2,Int64}(
  [1 0; 0 1],
  [0, 0]
)
"""
struct SpaceGroupElement{N,T<:Integer} <: GroupElement
    a::SMatrix{N,N,T}
    b::SVector{N, Rational{T}}
end

"""
    reduce(g::SpaceGroupElement{N,T}) where {N, T<:Integer}

Brings the translation vector of a space group element inside the standard unit cell.

# Arguments
- `g::SpaceGroupElement{N,T}`: The space group element to be reduced.

# Returns
- A new `SpaceGroupElement` with the same linear transformation matrix and a reduced translation vector.
"""
function reduce(g::SpaceGroupElement{N,T}) where{N, T<:Integer}
    return SpaceGroupElement{N,T}(g.a, g.b-floor.(g.b))
end

"""
    SpaceGroupElement{N,T}()

Constructor with no arguments, which creates the identity element of the space group.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

# Returns
- The identity element of the space group.
"""
SpaceGroupElement{N,T}() where {N,T<:Integer} =
    SpaceGroupElement(one(SMatrix{N,N,T}), zeros(SVector{N, Rational{T}}))

"""
    SpaceGroupElement{N,T}(t::SVector{N,T}) where {N, T<:Integer}

Constructor with a translation vector argument, which creates a pure translation element of the space group.

# Arguments
- `t::SVector{N,T}`: The translation vector.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

# Returns
- A space group element with the identity matrix as the linear transformation matrix and the given translation vector.    
"""
SpaceGroupElement(t::SVector{N,T})  where{N, T<:Integer} = 
    SpaceGroupElement(one(SMatrix{N,N,T}), Rational{T}.(t))

"""
    SpaceGroupElement{N,T}(m::SMatrix{N,N,T}) where {N, T<:Integer}

Constructor with a transformation matrix argument, which creates a pure linear transformation element of the space group.

# Arguments
- `m::SMatrix{N,N,T}`: The transformation matrix.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

# Returns
- A space group element with the given transformation matrix and zero translation vector.
"""
SpaceGroupElement(m::SMatrix{N,N,T})  where{N, T<:Integer} =
    SpaceGroupElement(m, zeros(SVector{N, Rational{T}}))

"""
    *(e1::SpaceGroupElement{N,T}, e2::SpaceGroupElement{N,T}) where {N,T<:Integer}

Composition of two space group elements.

# Arguments

- `e1::SpaceGroupElement{N,T}`: The first space group element.
- `e2::SpaceGroupElement{N,T}`: The second space group element.

# Returns
- The composition of the two space group elements.
"""
function *(e1::SpaceGroupElement{N,T}, e2::SpaceGroupElement{N,T}) where {N,T<:Integer}
    SpaceGroupElement(e1.a*e2.a, e1.a*e2.b+e1.b)
end

"""
    ∘(e1::SpaceGroupElement{N,T}, e2::SpaceGroupElement{N,T}) where {N,T<:Integer}

    Reduced composition of two space group elements. The translation vector of the result is brought 
    inside the standard unit cell.

    # Arguments
    - `e1::SpaceGroupElement{N,T}`: The first space group element.
    - `e2::SpaceGroupElement{N,T}`: The second space group element.

    # Returns
    - The reduced composition of the two space group elements.
"""
function ∘(e1::SpaceGroupElement{N,T}, e2::SpaceGroupElement{N,T}) where {N,T<:Integer}
    reduce(e1*e2)
end

"""
    SpaceGroupQuotient{N,T} = FiniteGroup{SpaceGroupElement{N,T}}

    The factor group of the space group with respect to the subgroup of pure translations.
    This group is isomorphic to the point group of the space group, but the representation 
    of its elements contains enough information to reconstruct the original space group. In 
    particular, the representation of symmorphic and non-symmorphic space groups is different.

    # Type Parameters
    - `N`: The dimension of the space.
    - `T<:Integer`: The type of the elements in the transformation matrix and translation vector.
"""
const SpaceGroupQuotient{N,T} = FiniteGroup{SpaceGroupElement{N,T}}

