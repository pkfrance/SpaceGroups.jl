import Base: length, eltype, iterate

"""
    GroupElement
    
An abstract type representing an element of a group.
"""
abstract type GroupElement end


"""
    FiniteGroup{E <: GroupElement}

Represents a finite group with elements of type `E`.

# Type Parameters
- `E <: GroupElement`: The type of the elements in the group, which must be a subtype of `GroupElement` and provide a composition operation `∘`.

# Fields
- `e::Set{E}`: The set of all group elements.

# Example
```julia
struct MyElement <: GroupElement
    # Define the structure of your group element here
end

function ∘(x::MyElement, y::MyElement)
    # Define the group operation here
end

# Define the generators of the group
gen = (MyElement(),)  # An iterable of group elements
group = FiniteGroup{MyElement}(gen)
```
"""
struct FiniteGroup{E <: GroupElement}
    e::Set{E}

    function (::Type{<:FiniteGroup})(gen)
        new{Base.eltype(gen)}(generate_(gen))
    end
end

"""
    FiniteGroup{E <: GroupElement}() where {E <: GroupElement}
    Default constructor creates the trivial finite group with only the identity element.
"""
FiniteGroup{E}() where {E <: GroupElement} = FiniteGroup{E}(Set{E}())

function generate_(gen)
    identity=eltype(gen)()
    e = Set((identity,)) 
    last = e
    while true
        found = setdiff((x ∘ y for x in gen, y in last), e)
        if isempty(found) break end
        union!(e, found)
        last = found
    end
    e
end

# Iteration
iterate(fg::FiniteGroup) = iterate(fg.e)
iterate(fg::FiniteGroup, i) = iterate(fg.e, i)

# Length and element type
eltype(fg::FiniteGroup) = eltype(fg.e)
length(fg::FiniteGroup) = length(fg.e)
