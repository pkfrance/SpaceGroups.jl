

"""
    AffinePhase{N, T<:Integer}

Represents the phase of a plane wave x ↦ exp(2π i (k ⋅ x + ϕ)).

The wave is characterized by a wave vector `k` and a phase offset `ϕ`.

# Type Parameters
- `N`: The dimensionality of the space.
- `T<:Integer`: The numeric type for the elements of the wave vector.

# Fields
- `k::SVector{N, T}`: The wave vector.
- `ϕ::Rational{T}`: The phase offset, always maintained within the interval `[0, 1)`.

# Constructors
- `AffinePhase(k::SVector{N,T}, ϕ::Rational{T})`: Creates a new `AffinePhase` instance. 
The phase `ϕ` is automatically normalized to its fractional part.
"""
struct AffinePhase{N,T<:Integer}
    k::SVector{N,T}
    ϕ::Rational{T}
    AffinePhase{N,T}(k::SVector{N,T}, ϕ::Rational{T}) where {N,T<:Integer} = new(k, mod(ϕ, 1))
end

AffinePhase(k::SVector{N,T}, ϕ::Rational{T}) where {N,T<:Integer} =
    AffinePhase{N,T}(k, ϕ)


function Base.show(io::IO, ::MIME"text/plain", x::AffinePhase)
    print(io, "AffinePhase(k = $(x.k), ϕ = $(x.ϕ))")
end

function Base.show(io::IO, x::AffinePhase)
    print(io, "AP($(x.k), $(x.ϕ))")
end

"""
    *(g::SpaceGroupElement, ap::AffinePhase) -> AffinePhase

Applies a symmetry operation `g` to an `AffinePhase` `ap`.

The action of a space group element `g = (a, b)` on an `AffinePhase` with wave vector `k` and phase `ϕ`.
This results in a new `AffinePhase` with a transformed wave vector `k' = aᵀk` and a new phase `ϕ' = ϕ + b ⋅ k`.

# Arguments
- `g::SpaceGroupElement{N,T}`: The symmetry operation to apply.
- `ap::AffinePhase{N,T}`: The affine phase to be transformed.

# Returns
- `AffinePhase{N,T}`: A new `AffinePhase` instance after the transformation.
"""
function *(g::SpaceGroupElement{N,T}, ap::AffinePhase{N,T}) where {N,T<:Integer}
    # Apply the symmetry operation to the phase
    AffinePhase(transpose(g.a) * ap.k, ap.ϕ + g.b ⋅ ap.k)
end

"""
    ComplexOrbit{N, T<:Integer}

Represents an orbit of plane waves where the phase is unconstrained.

In a `ComplexOrbit`, for any wave vector `k` in the orbit, its antipode `-k` is not present.

# Type Parameters
- `N`: The dimensionality of the space.
- `T<:Integer`: The numeric type for the elements of the affine phase.

# Fields
- `aps::Vector{AffinePhase{N,T}}`: A vector of `AffinePhase` objects that constitute the orbit.
"""
struct ComplexOrbit{N,T<:Integer}
    aps::Vector{AffinePhase{N,T}}
end

"""
    RealOrbit{N, T<:Integer}

Represents an orbit of plane waves where the phase is fixed modulo π.

In a `RealOrbit`, for every wave vector `k` in the orbit, its antipode `-k` is also present. 
The phase is specifically determined to ensure that the superposition of two antipodal waves 
results in a real-valued function.

# Type Parameters
- `N`: The dimensionality of the space.
- `T<:Integer`: The numeric type for the elements of the affine phase.

# Fields
- `aps::Vector{AffinePhase{N,T}}`: A vector of `AffinePhase` objects forming the orbit.
"""
struct RealOrbit{N,T<:Integer}
    aps::Vector{AffinePhase{N,T}}
end

"""
    ExtinctOrbit{N, T<:Integer}

Represents an orbit of wave vectors that corresponds to a systematic extinction.

Systematic extinctions arise in the case of non-symmorphic space groups when a symmetry 
operation preserves the direction of a wave vector but alters the phase of the corresponding 
plane wave, leading to destructive interference.

# Type Parameters
- `N`: The dimensionality of the space.
- `T<:Integer`: The numeric type for the elements of the affine phase.

# Fields
- `k::Vector{SVector{N,T}}`: A vector of wave vectors `k` that are subject to extinction.
"""
struct ExtinctOrbit{N,T<:Integer}
    k::Vector{SVector{N,T}}
end


"""
    FormalOrbit{N, T<:Integer}

A type alias representing the union of `ComplexOrbit`, `RealOrbit`, and `ExtinctOrbit`.

This union type encompasses all possible outcomes of generating an orbit from a wave vector and a 
space group, including orbits that may be physically unobservable due to extinction.
"""
const FormalOrbit{N,T<:Integer} = Union{ComplexOrbit{N,T},RealOrbit{N,T},ExtinctOrbit{N,T}}


"""
    PhysicalOrbit{N, T<:Integer}

A type alias for the union of `ComplexOrbit` and `RealOrbit`.

This type represents orbits corresponding to actual Bragg peaks.
"""
const PhysicalOrbit{N,T<:Integer} = Union{ComplexOrbit{N,T},RealOrbit{N,T}}


"""
    make_orbit(k::SVector{N,T}, G::SpaceGroupQuotient{N,T}) -> FormalOrbit{N,T}

Generates an orbit of affine phases from an initial wave vector `k` for a given 
space group (quotient) `G`.

The function determines the type of orbit (`Complex`, `Real`, or `Extinct`) by 
analyzing the action of the symmetry operations in `G` on the initial wave vector.

# Arguments
- `k::SVector{N,T}`: The starting wave vector.
- `G::SpaceGroupQuotient{N,T}`: The quotient of the space group.

# Returns
- `FormalOrbit{N,T}`: The generated orbit, which can be a `ComplexOrbit`, `RealOrbit`, or `ExtinctOrbit`.

# Logic
1.  **Extinction Check**: If any symmetry operation maps `k` to itself but shifts its phase, the orbit is `ExtinctOrbit`.
2.  **Real Orbit Check**: If the orbit contains the antipode `-k` for any `k` in the orbit, it is a `RealOrbit`. In this case, 
    for the corresponding affine phases, ϕ is shifted to have opposite sign for the antipodes to insure that their 
    superposition is a real function.  
3.  **Complex Orbit**: If neither of the above conditions is met, the orbit is `ComplexOrbit`.

# Example
```julia-repl
julia> g=@SGE([-1 0; 0 1], [0//1, 1//2]);

julia> p1g1=SpaceGroupQuotient([g])
SpaceGroupQuotient (dimension 2, order 2)

julia> make_orbit([1, -1], p1g1)
ComplexOrbit with 2 elements

julia> make_orbit([1, 0], p1g1)
RealOrbit with 2 elements

julia> make_orbit([0, 1], p1g1)
ExtinctOrbit with 2 elements
```
"""
function make_orbit(k::SVector{N,T}, G::SpaceGroupQuotient{N,T})::FormalOrbit{N,T} where {N,T<:Integer}
    d = Dict{SVector{N,T},Rational{T}}()
    isextinct = false
    isreal = false
    ap0 = AffinePhase(k, zero(Rational{T}))
    for e in G
        ap = e * ap0
        d[ap.k] = ap.ϕ
        if ap.k == k && !iszero(ap.ϕ)
            # If the symmetry operation produces the same wave-vector 
            # with a different phase, the orbit is extinct
            isextinct = true
        end
        if ap.k == -k
            # If the orbit contains the antipode, it is real
            isreal = true
        end
    end
    if isextinct
        ks = [k for k in keys(d)]
        if isreal # The antipodes are already there
            return ExtinctOrbit(ks)
        else
            ksa = [-k for k in keys(d)] # Add antipodes
            return ExtinctOrbit(vcat(ks, ksa))
        end
    end
    if isreal
        # The phase of each partial wave is the average 
        # of the phases of the antipodal waves:
        return RealOrbit([AffinePhase(k, (d[k] - d[-k]) // 2) for k in keys(d)])
    else
        # If the orbit is not extinct and contains no antipodes,
        # the phase of each partial wave is already known up to a global phase:
        return ComplexOrbit([AffinePhase(k, d[k]) for k in keys(d)])
    end
end

# A convenience version of `make_orbit` that accepts an arbitrary vector type for `k`
# instead of a stativ vector only.
function make_orbit(k::AbstractVector{T}, G::SpaceGroupQuotient{N,T}) where {N,T<:Integer}
    length(k) == N || throw(ArgumentError("Expected vector of length $N, got length $(length(k))"))
    make_orbit(SVector{N,T}(k), G)
end


function Base.show(io::IO, ::MIME"text/plain", x::ComplexOrbit)
    print(io, "ComplexOrbit with $(length(x.aps)) elements")
end

function Base.show(io::IO, x::ComplexOrbit)
    print(io, "ComplexOrbit with $(length(x.aps)) elements")
end

function Base.show(io::IO, ::MIME"text/plain", x::RealOrbit)
    print(io, "RealOrbit with $(length(x.aps)) elements")
end

function Base.show(io::IO, x::RealOrbit)
    print(io, "RealOrbit with $(length(x.aps)) elements")
end

function Base.show(io::IO, ::MIME"text/plain", x::ExtinctOrbit)
    print(io, "ExtinctOrbit with $(length(x.k)) elements")
end

function Base.show(io::IO, x::ExtinctOrbit)
    print(io, "ExtinctOrbit with $(length(x.k)) elements")
end