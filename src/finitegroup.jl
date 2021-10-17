import Base: length, eltype

abstract type GroupElement end

struct FiniteGroup{E <: GroupElement}
    e::Set{E}

    function (::Type{<:FiniteGroup})(gen)
        new{Base.eltype(gen)}(generate_(gen))
    end
end

function generate_(gen)
    e = Set(gen)
    last = e
    while true
        found = setdiff((x * y for x in gen, y in last), e)
        if isempty(found) break end
        union!(e, found)
        last = found
    end
    e
end

# Iteration
iterate(fg::FiniteGroup) = iterate(fg.e)
iterate(fg::FiniteGroup, i) = iterate(fg.e, i)
eltype(fg::FiniteGroup) = eltype(fg.e)
length(fg::FiniteGroup) = length(fg.e)