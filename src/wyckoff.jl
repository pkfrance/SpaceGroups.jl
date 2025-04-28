import Base.*

"""
    WyckoffPosition{N,T<:Integer}

A Wyckoff position in N dimensions.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

# Fields
- `anchor::StaticArrays.SVector{N, Rational{T}}`: The anchor point of the Wyckoff position.
- `directions::StaticArrays.SMatrix{N,M,T}`: The directions of the Wyckoff position.

# Constructors
- `WyckoffPosition{N,T}(anchor::SVector{N,Rational{T}}, directions::SMatrix{N,M,T})`: A Wyckoff position with 
free parameters. The constructore checks that the directions are linearly independent. 
- `WyckoffPosition{N,T}()`: A generic Wyckoff position, the number of free parameters equals the dimension of the space. 
- `WyckoffPosition{N,T}(anchor::SVector{N,Rational{T}})`: A Wyckoff position with zero free parameters (the most 
special type of Wyckoff position).
"""
struct WyckoffPosition{N, T<:Integer}
    anchor::SVector{N, Rational{T}}
    directions::SMatrix{N, M, T} where M
    function WyckoffPosition(anchor::SVector{N, Rational{T}}, directions::SMatrix{N, M, T}) where {N, M, T<:Integer}
        if rank(directions) !=M
            throw(ArgumentError("The directions must be linearly independent."))
        end
        new{N, T}(anchor, directions)
    end
end

# Default constructor creates a generic Wyckoff position
function WyckoffPosition{N, T}() where {N, T<:Integer}
    WyckoffPosition(zero(SVector{N, T}), SMatrix{N,N,T}(I))
end


# Constructor in the case of zero free parameters of the Wyckoff position
"""
    WyckoffPosition{N,T}(anchor::SVector{N,Rational{T}}) where {N, T<:Integer}
Create a Wyckoff position with zero free parameters.
# Example
```julia
julia> using StaticArrays
julia> w = WyckoffPosition{2,Int}(SA[1//2, 1//2])
WyckoffPosition{2, Int64}(Rational{Int64}[1//2, 1//2], 2×0 SMatrix{2, 0, Int64, 0} with indices SOneTo(2)×SOneTo(0))

julia> w.anchor
2-element SVector{2, Rational{Int64}} with indices SOneTo(2):
 1//2
 1//2

julia> w.directions
2×0 SMatrix{2, 0, Int64, 0} with indices SOneTo(2)×SOneTo(0)
```
"""
function WyckoffPosition{N, T}(anchor::SVector{N, Rational{T}}) where {N, T<:Integer}
    WyckoffPosition(anchor, SMatrix{N,0,T}())
end

"""
    *(e::SpaceGroupElement{N,T}, w::WyckoffPosition{N, T}) -> WyckoffPosition{N, T}

Action of a space group element on a Wyckoff position.

# Example
```julia
julia> using StaticArrays

julia> e = SpaceGroupElement{2,Int}(SMatrix{2,2,Int}([1 0; 0 1]))
SpaceGroupElement{2,Int64}(
  [1 0; 0 1],
  [0, 0]
)

julia> w = WyckoffPosition{2,Int}(SVector{2,Int}([1, 1]))
WyckoffPosition{2,Int64}(
  [1, 1],
  Int64[]
)

julia> e*w
WyckoffPosition{2,Int64}(
  [1, 1],
  Int64[]
)
```
"""
function *(e::SpaceGroupElement{N,T}, w::WyckoffPosition{N, T}) where {N, T<:Integer}
    anchor=e.a*w.anchor+e.b
    directions=e.a*w.directions
    WyckoffPosition{N,T}(anchor, directions)
end

"""
    normalize(w::WyckoffPosition{N, T}) -> Tuple{WyckoffPosition{N, T}, StaticArrays.SVector{N, T}}

Normalize a Wyckoff position to the standard unit cell.

# Example
```julia
julia> using StaticArrays

julia> w = WyckoffPosition{2,Int}(SVector{2,Int}([1, 1]))
WyckoffPosition{2,Int64}(
  [1, 1],
  Int64[]
)

julia> normalize(w)
(WyckoffPosition{2,Int64}(
  [1, 1],
  Int64[]
), [1, 1])
```
"""
function normalize(w::WyckoffPosition{N, T})::Tuple{WyckoffPosition{N, T}, SVector{N, T}} where {N, T<:Integer}
    t=T.(floor.(w.anchor))
    (WyckoffPosition(w.anchor-t, w.directions), t)
end

"""
    stabilizer_quotient(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T}) where {N, T<:Integer}
Compute the quotient of the stabilizer group of a Wyckoff position (with respect to translations).
# Arguments
- `w::WyckoffPosition{N, T}`: The Wyckoff position.
- `G::SpaceGroupQuotient{N, T}`: The space group quotient.
# Returns
- The stabilizer subgroup of G for the Wyckoff position w.
"""
function stabilizer_quotient(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T})::SpaceGroupQuotient{N, T} where {N, T<:Integer}
    # Compute the quotient of the stabilizer group of a Wyckoff position with respect to translations
    s=Set{SpaceGroupElement{N, T}}()
    for g in G
        if normalize(g*w) == normalize(w) # Select element acting on w by a lattice translation
            push!(s, g)
        end
    end
    return SpaceGroupQuotient(s)
end

"""
    is_valid(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T})::Bool where {N, T<:Integer}
Check if `w` is a valid special Wyckoff position for the space group quotient G. Namely, check 
if the stabilizer group of `w` does not conserve more directions than the number of free parameters of `w`.
# Arguments
- `w::WyckoffPosition{N, T}`: The Wyckoff position.
- `G::SpaceGroupQuotient{N, T}`: The space group quotient.
# Returns
- `true` if `w` is a valid special Wyckoff position for the space group quotient `G`, `false` otherwise.
"""
function is_valid(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T})::Bool where {N, T<:Integer}
    s=[g.a-I for g in stabilizer_quotient(w, G)]
    # Check the dimension of the common kernel of all elements of `s`
    kernel_dim=N-rank(vcat(s...))
    if kernel_dim != size(w.directions, 2)
        return false
    end
    true
end
