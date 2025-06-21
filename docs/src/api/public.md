# Public API

```@meta
CurrentModule = SpaceGroups
```

## Algebra
```@docs
SpaceGroupElement
SpaceGroupQuotient
@SGE
Base.:*(::SpaceGroupElement{N,T}, ::SpaceGroupElement{N,T}) where {N,T<:Integer}
∘(::SpaceGroupElement{N,T}, ::SpaceGroupElement{N,T}) where {N,T<:Integer}
reduce
```

## Reciprocal space 
```@docs
AffinePhase
Base.:*(::SpaceGroupElement{N,T}, ::AffinePhase{N, T}) where {N,T<:Integer}
ComplexOrbit
RealOrbit
ExtinctOrbit
FormalOrbit
PhysicalOrbit
make_orbit
```

## Wyckoff positions
```@docs
WyckoffPosition
@WP
Base.:*(::SpaceGroupElement{N,T}, ::WyckoffPosition{N, T}) where {N,T<:Integer}
stabilizer_quotient
is_valid_wyckoff
normalize
```