"""
    SpaceGroupQuotient{N,T} = FiniteGroup{SpaceGroupElement{N,T}}

The factor group of the space group with respect to the subgroup of pure translations.
This group is isomorphic to the point group of the space group, but the representation 
of its elements contains enough information to reconstruct the original space group. In 
particular, the representation of symmorphic and non-symmorphic space groups is different.

# Type Parameters
- `N`: The dimension of the space.
- `T<:Integer`: The type of the elements in the transformation matrix and translation vector.

# Constructors
- SpaceGroupQuotient{N,T}(): Construct the trivial space group (P1) of dimension N 
- SpaceGroupQuotient{N,T}(gen::AbstractVector{SpaceGroupElement{N,T}}): Constructs the space 
  group using the generating set `gen` 

# Examples
```julia-repl
julia> SpaceGroupQuotient{2, Int}()
SpaceGroupQuotient (dimension 2, order 1)

julia> g1=@SGE([-1 0; 0 -1])
SpaceGroupElement(
  a = [-1 0; 0 -1],
  b = [0//1, 0//1]
)

julia> g2=@SGE([-1 0; 0 1], [1//2, 0//1])
SpaceGroupElement(
  a = [-1 0; 0 1],
  b = [1//2, 0//1]
)

julia> p2mg=SpaceGroupQuotient([g1, g2])
SpaceGroupQuotient (dimension 2, order 4)
```
"""
struct SpaceGroupQuotient{N,T<:Integer} 
    group::FiniteGroup{SpaceGroupElement{N,T}}
end

SpaceGroupQuotient{N, T}() where {N,T<:Integer} = 
    SpaceGroupQuotient{N, T}(FiniteGroup{SpaceGroupElement{N,T}}())

function SpaceGroupQuotient(gen::AbstractVector{SpaceGroupElement{N,T}}) where {N,T<:Integer}
    group=FiniteGroup{SpaceGroupElement{N,T}}(gen)
    SpaceGroupQuotient{N, T}(group)
end

@delegate SpaceGroupQuotient.group (Base.length, Base.iterate, Base.eltype)

function Base.show(io::IO, ::MIME"text/plain", G::SpaceGroupQuotient{N,T}) where {N,T}
    print(io, "SpaceGroupQuotient (dimension $N, order $(length(G)))")
end
