import Base: length, eltype, iterate

"""
    GroupElement
    
An abstract type representing an element of a group.
"""
abstract type GroupElement end

function generate_(gen::AbstractVector{E})::Set{E} where E<:GroupElement
    identity=E()
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


struct FiniteGroup{E <: GroupElement}
    d::Dict{E, Int}

    function (::Type{<:FiniteGroup})(gen::AbstractVector{E}) where E <: GroupElement
        e=generate_(gen)
        identity=E()
        # By convention, the index of the neutral element is 1
        d=Dict{E, Int}(identity=>1) 
        i=1
        for k in e
            if k!=identity
                i+=1
                d[k]=i
            end
        end
        new{E}(d)
    end
end

"""
    FiniteGroup{E <: GroupElement}() where {E <: GroupElement}
    Default constructor creates the trivial finite group with only the identity element.
"""
FiniteGroup{E}() where {E <: GroupElement} = FiniteGroup{E}(E[])


# Iteration
iterate(fg::FiniteGroup) = iterate(keys(fg.d))
iterate(fg::FiniteGroup, i) = iterate(keys(fg.d), i)

# Length and element type
eltype(fg::FiniteGroup) = eltype(keys(fg.d))
length(fg::FiniteGroup) = length(fg.d)
