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
- `WyckoffPosition{N,T}(anchor::SVector{N,Rational{T}})`: A Wyckoff position with zero directions.
"""
struct WyckoffPosition{N, T<:Integer}
    anchor::SVector{N, Rational{T}}
    directions::SMatrix{N, M, T} where M
end

# Constructor in the case of zero parameters
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