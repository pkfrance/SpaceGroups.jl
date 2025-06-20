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
- `WyckoffPosition(anchor::SVector{N,Rational{T}}, directions::SMatrix{N,M,T})`: A Wyckoff 
  position with free parameters. The constructore checks that the directions are 
  linearly independent. 
- `WyckoffPosition{N,T}()`: A general Wyckoff position, the number of free parameters 
  equals the dimension of the space. 
- `WyckoffPosition(anchor::SVector{N,Rational{T}})`: A Wyckoff position with zero 
  free parameters (the most special type of Wyckoff position).
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

# Default constructor creates a general Wyckoff position
function WyckoffPosition{N, T}() where {N, T<:Integer}
    WyckoffPosition(zero(SVector{N, Rational{T}}), SMatrix{N,N,T}(I))
end

# Convenience default constructor
WyckoffPosition{N}() where N = WyckoffPosition{N, Int}()


# Constructor in the case of zero free parameters of the Wyckoff position
"""
    WyckoffPosition{N,T}(anchor::SVector{N,Rational{T}}) where {N, T<:Integer}
Create a Wyckoff position with zero free parameters.
# Example
```julia-repl
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
function WyckoffPosition(anchor::SVector{N, Rational{T}}) where {N, T<:Integer}
    WyckoffPosition(anchor, SMatrix{N,0,T}())
end

"""
    @WP(anchor_expr, dirs_expr=nothing)
Construct a Wyckoff position using a convenient macro interface.

- `anchor_expr`: An expression that evaluates to a vector of rationals, specifying the anchor point.
- `dirs_expr`: (Optional) An expression that evaluates to a matrix of integers, specifying the directions. 
  If omitted or `nothing`, the Wyckoff position will have zero free parameters.

# Examples
```julia-repl
julia> @WP([0//1, 1//2, 1//2], [1; 1; 1;;])
WyckoffPosition(anchor=[0//1, 1//2, 1//2]), directions=[1; 1; 1;;])

julia> @WP([0//1, 1//2, 1//2])
WyckoffPosition(anchor=[0//1, 1//2, 1//2]), no parameters)
```
"""
macro WP(anchor_expr, dirs_expr=nothing)
    quote
        # Evaluate anchor vector
        anchor_val = $(esc(anchor_expr))
        if !(anchor_val isa AbstractVector{<:Rational})
            throw(ArgumentError("Anchor must be a vector of Rationals."))
        end

        N = length(anchor_val)
        T = promote_type(typeof(numerator(anchor_val[1])), typeof(denominator(anchor_val[1])))

        # If only anchor is given, use empty direction matrix
        if $(dirs_expr) === nothing
            directions = SMatrix{N, 0, T}(reshape(T[], N, 0))
            WyckoffPosition(SVector{N, Rational{T}}(anchor_val), directions)
        else
            dirs_val = $(esc(dirs_expr))
            if !(dirs_val isa AbstractMatrix{<:Integer})
                throw(ArgumentError("Directions must be a matrix of Integers."))
            end

            # Check dimension consistency
            size(dirs_val, 1) == N || throw(DimensionMismatch("Anchor and direction matrix row counts must match."))

            # Promote integer type
            T2 = promote_type(T, eltype(dirs_val))
            anchor_promoted = SVector{N, Rational{T2}}(Rational{T2}.(anchor_val))
            M = size(dirs_val, 2)
            dirs_promoted = SMatrix{N, M, T2}(T2.(dirs_val))

            WyckoffPosition(anchor_promoted, dirs_promoted)
        end
    end
end


"""
    *(e::SpaceGroupElement{N,T}, w::WyckoffPosition{N, T}) -> WyckoffPosition{N, T}

Action of a space group element on a Wyckoff position.

# Example
```julia-repl
julia> g=@SGE([0 1; -1 0])
SpaceGroupElement(
  a = [0 1; -1 0],
  b = [0//1, 0//1]
)

julia> w=@WP([1//2, 0//1])
WyckoffPosition(anchor=[1//2, 0//1]), no parameters)

julia> g*w
WyckoffPosition(anchor=[0//1, -1//2]), no parameters)
```
"""
function *(e::SpaceGroupElement{N,T}, w::WyckoffPosition{N, T}) where {N, T<:Integer}
    anchor=e.a*w.anchor+e.b
    directions=e.a*w.directions
    WyckoffPosition(anchor, directions)
end

"""
    normalize(w::WyckoffPosition{N, T}) -> Tuple{WyckoffPosition{N, T}, StaticArrays.SVector{N, T}}

Normalize a Wyckoff position to the standard unit cell.

# Example
```julia-repl
julia> w=@WP([1//1, 1//1])
WyckoffPosition(anchor=[1//1, 1//1]), no parameters)

julia> normalize(w)
(WP{2,0}, [1, 1])
```
"""
function normalize(w::WyckoffPosition{N, T})::Tuple{WyckoffPosition{N, T}, SVector{N, T}} where {N, T<:Integer}
    t=T.(floor.(w.anchor))
    (WyckoffPosition(w.anchor-t, w.directions), t)
end

"""
    stabilizer_quotient(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T}) where {N, T<:Integer}
Compute the quotient of the stabilizer group (the site symmetry group) of a Wyckoff position 
(with respect to translations).
# Arguments
- `w::WyckoffPosition{N, T}`: The Wyckoff position.
- `G::SpaceGroupQuotient{N, T}`: The space group quotient.
# Returns
- The stabilizer subgroup of G for the Wyckoff position w.
# Example
```julia-repl
julia> g1=@SGE([-1 0; 0 -1]);

julia> g2=@SGE([-1 0; 0 1], [1//2, 0//1]);

julia> p2mg=SpaceGroupQuotient([g1, g2]);

julia> w1=@WP([1//4, 0], [0; 1;;]);

julia> G=stabilizer_quotient(w1, p2mg)
SpaceGroupQuotient (dimension 2, order 2)

julia> G.e
Set{SpaceGroupElement{2, Int64}} with 2 elements:
  SGE([1 0; 0 1], [0//1, 0//1])
  SGE([-1 0; 0 1], [1//2, 0//1])
```
"""
function stabilizer_quotient(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T})::SpaceGroupQuotient{N, T} where {N, T<:Integer}
    # Compute the quotient of the stabilizer group of a Wyckoff position with respect to translations
    s=Set{SpaceGroupElement{N, T}}()
    w0, _ = normalize(w) 
    for g in G
        u, _= normalize(g*w0)
        if u == w0 # Check if the normalized Wyckoff position is invariant under the group element
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
# Example
```julia-repl
julia> g1=@SGE([-1 0; 0 -1]);

julia> g2=@SGE([-1 0; 0 1], [1//2, 0//1]);

julia> p2mg=SpaceGroupQuotient([g1, g2]);)

julia> w1=@WP([1//4, 0], [0; 1;;])
WyckoffPosition(anchor=[1//4, 0//1]), directions=[0; 1;;])

julia> is_valid_wyckoff(w1, p2mg)
true

julia> w2=@WP([1//3, 0], [0; 1;;])
WyckoffPosition(anchor=[1//3, 0//1]), directions=[0; 1;;])

julia> is_valid_wyckoff(w2, p2mg)
false
```
"""
function is_valid_wyckoff(w::WyckoffPosition{N, T}, G::SpaceGroupQuotient{N, T})::Bool where {N, T<:Integer}
    s=[g.a-I for g in stabilizer_quotient(w, G)]
    # Check the dimension of the common kernel of all elements of `s`
    kernel_dim=N-rank(vcat(s...))
    if kernel_dim != size(w.directions, 2)
        return false
    end
    true
end

import Base: show, summary

# Long form: used in REPL and when printed directly
function Base.show(io::IO, ::MIME"text/plain", wp::WyckoffPosition{N, T}) where {N, T}
    M = size(wp.directions, 2)

    if M == N
        print(io, "WyckoffPosition(dim=$N, general)")
    else
        print(io, "WyckoffPosition(")
        print(io, "anchor=[")
        for (i, val) in enumerate(wp.anchor)
            i > 1 && print(io, ", ")
            print(io, val)  # print each Rational cleanly
        end
        print(io, "])")    
        if M == 0
            print(io, ", no parameters")
        else
            print(io, ", directions=", wp.directions)
        end
        print(io, ")")
    end
end

# Short form: used in collections
function Base.show(io::IO, wp::WyckoffPosition{N, T}) where {N, T}
    M = size(wp.directions, 2)
    print(io, "WP{$N,$M}")
end
