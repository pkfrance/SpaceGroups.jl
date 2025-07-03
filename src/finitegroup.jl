import Base: length, eltype, iterate

"""
    GroupElement
    
An abstract type representing an element of a group.
"""
abstract type GroupElement end

# Fallback methods to enforce interface implementation
identity(::Type{E}) where E<:GroupElement = error("Define identity(::Type{$E})")

function ∘(a::E, b::E) where E<:GroupElement
    msg = """
    The group operation '∘' not implemented for group elements of type `$E`.
    To use this type with the FiniteGroup structure, define the method:
    `∘(a::$E, b::$E) = ...`
    """
    error(msg)
end

function generate_(gen::AbstractVector{E})::Set{E} where E<:GroupElement
    s = Set((identity(E),))
    last = s
    while true
        found = setdiff((x ∘ y for x in gen, y in last), s)
        if isempty(found)
            break
        end
        union!(s, found)
        last = found
    end
    s
end


struct OperationCache
    tabmul::Matrix{Int}
    tabinv::Vector{Int}
end

struct ConjugacyCache
    conjind::Vector{Int}
    conjcomp::Vector{Vector{Int}}
end

mutable struct FiniteGroupCache
    oc::Union{Nothing,OperationCache}
    cc::Union{Nothing, ConjugacyCache}

    FiniteGroupCache() = new(nothing, nothing)
end


"""
    FiniteGroup{E <: GroupElement}

Represents a finite group with elements of type `E`.

# Type Parameters
- `E <: GroupElement`: The type of the elements in the group, which must be a 
  subtype of `GroupElement` and provide the group operation `∘`, which must be associative, as well as the 
  function `identity(::Type{E})` returning the neutral element of the group. 
  The type `E` should be preferrably immutable, as it will be the type of keys 
  in a dictionary. The equality of objects of type `E` should be equivalent to 
  equality of the corresponding group elements. 

# Fields
- `d::Dict{E, Int}`: The keys of `d` are group elements and the values are 
  arbitrarily assigned indices (with the convention that the neutral element
  corresponds to the index 1). Reciprocal to `e`: `d[e[i]]=i`
- `e::Vector{E}`: The vector of group elements. Reciprocal to `d`: `e[d[g]]=g` 
- `_cache::FiniteGroupCache{E}`: This structure holds all data derived from `d`,
  and lazily inititialized when the corresponding functionality is requested
- `_lock::ReentrantLock`: Lock for thread-safe lazy initialization

# Example
```julia
struct MyElement <: GroupElement
    # Define the structure of your group element here
end

function ∘(x::MyElement, y::MyElement)
    # Define the group operation here
end

# Define the generators of the group
gen = [MyElement()]  # An `AbstractVector{MyElement}`
group = FiniteGroup{MyElement}(gen)
```
"""
struct FiniteGroup{E<:GroupElement}
    d::Dict{E,Int}
    e::Vector{E}
    _cache::FiniteGroupCache
    _lock::ReentrantLock

    function (::Type{<:FiniteGroup})(gen::AbstractVector{E}) where E<:GroupElement
        s = generate_(gen)
        id = identity(E)
        e = Vector{E}(undef, length(s))
        e[1] = id # By convention, the index of the neutral element is 1
        e[2:end] .= (g for g in s if g != id)
        d = Dict{E,Int}(g => i for (i, g) in enumerate(e))
        i = 1
        new{E}(d, e, FiniteGroupCache(), ReentrantLock())
    end
end

"""
    FiniteGroup{E <: GroupElement}() where {E <: GroupElement}
    Default constructor creates the trivial finite group with only the identity element.
"""
FiniteGroup{E}() where {E<:GroupElement} = FiniteGroup{E}(E[])





function _ensure_cache(G::FiniteGroup, field::Symbol, fieldtype::DataType)
    isnothing(getfield(G._cache, field)) || return nothing
    lock(G._lock) do
        isnothing(getfield(G._cache, field)) || return nothing
        setfield!(G._cache, field, fieldtype(G))
    end
    nothing
end



function OperationCache(G::FiniteGroup)
    n = length(G.e)
    tabmul = Matrix{Int}(undef, n, n)
    tabinv = Vector{Int}(undef, n)
    for i in 1:n
        for j in 1:n
            k = G.d[G.e[i]∘G.e[j]]
            tabmul[i, j] = k
            if k == 1
                tabinv[i] = j
            end
        end
    end
    OperationCache(tabmul, tabinv)
end



function _tabmul(G::FiniteGroup)::Matrix{Int}
    _ensure_cache(G, :oc, OperationCache)
    G._cache.oc.tabmul
end

function _tabinv(G::FiniteGroup)::Vector{Int}
    _ensure_cache(G, :oc, OperationCache)
    G._cache.oc.tabinv
end

function inv(g::E, G::FiniteGroup{E})::E where E<:GroupElement
    tabinv = _tabinv(G)
    G.e[tabinv[G.d[g]]]
end



function ConjugacyCache(G::FiniteGroup)
    tabmul = _tabmul(G)
    tabinv = _tabinv(G)
    conjind = zeros(Int, length(G))
    nc = 0
    for i in eachindex(conjind)
        if conjind[i] != 0
            continue
        end
        nc += 1 # New class found
        conjind[i] = nc
        for j in eachindex(G.e)
            jinv = tabinv[j]
            k = tabmul[jinv, tabmul[i, j]]
            conjind[k] = nc
        end
    end
    conjcomp = [Int[] for n in 1:nc]
    for i in eachindex(conjind)
        push!(conjcomp[conjind[i]], i)
    end
    ConjugacyCache(conjind, conjcomp)
end

function _conjind(G::FiniteGroup)::Vector{Int}
    _ensure_cache(G, :cc, ConjugacyCache)
    G._cache.cc.conjind
end

function _conjcomp(G::FiniteGroup)::Vector{Vector{Int}}
    _ensure_cache(G, :cc, ConjugacyCache)
    G._cache.cc.conjcomp
end

function conjugacy_class(g::E, G::FiniteGroup{E})::Vector{E} where E<:GroupElement
    conjind=_conjind(G)
    conjcomp=_conjcomp(G)
    ind=conjind[G.d[g]] # Class index
    [G.e[i] for i in conjcomp[ind]]
end

# Iteration
iterate(fg::FiniteGroup) = iterate(keys(fg.d))
iterate(fg::FiniteGroup, i) = iterate(keys(fg.d), i)

# Length and element type
eltype(fg::FiniteGroup) = eltype(keys(fg.d))
length(fg::FiniteGroup) = length(fg.d)
