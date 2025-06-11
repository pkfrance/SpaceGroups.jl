

"""
    AffinePhase{N, T}(k::SVector{N,T}, ϕ::Rational{T}) where {N, T<:Integer}

    Represents a phase of the planary wave with a wave vector `k` and a phase `2π⋅ϕ`.

    # Type Parameters

    - `N`: The dimension of the space.
    - `T<:Integer`: The type of the elements in the wave vector.

    # Fields

    - `k::StaticArrays.SVector{N,T}`: The wave vector.
    - `ϕ::Rational{T}`: The phase of the wave.

    # Constructors

    - `AffinePhase{N, T}(k::SVector{N,T}, ϕ::Rational{T})`: Create a new `AffinePhase` object 
    with the phase equal to the fractional part of ϕ.
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
    *(g::SpaceGroupElement, ap::AffinePhase{N,T}) where {N,T<:Integer}

Applies a symmetry operation represented by a `SpaceGroupElement` to an `AffinePhase`.

# Arguments
- `g::SpaceGroupElement`: The symmetry operation to be applied.
- `ap::AffinePhase{N,T}`: The affine phase to which the symmetry operation is applied. 
  Here, `N` represents the dimension of the space, and `T` is an integer type parameter.

# Returns
- `AffinePhase`: A new `AffinePhase` object resulting from applying the symmetry operation. 
  The phase vector `k` is transformed by the transpose of the matrix `a` from `g`, 
  and the phase offset `ϕ` is updated by adding the dot product of `k` and the translation vector `b` from `g`.
"""
function *(g::SpaceGroupElement{N,T}, ap::AffinePhase{N,T}) where {N,T<:Integer}
    # Apply the symmetry operation to the phase
    AffinePhase(transpose(g.a) * ap.k, ap.ϕ + g.b ⋅ ap.k)
end

struct ComplexOrbit{N,T<:Integer}
    aps::Vector{AffinePhase{N,T}}
end

struct RealOrbit{N,T<:Integer}
    aps::Vector{AffinePhase{N,T}}
end

struct ExtinctOrbit{N,T<:Integer}
    k::Vector{SVector{N,T}}
end

const FormalOrbit{N,T<:Integer} = Union{ComplexOrbit{N,T},RealOrbit{N,T},ExtinctOrbit{N,T}}

const PhysicalOrbit{N,T<:Integer} = Union{ComplexOrbit{N,T},RealOrbit{N,T}}

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