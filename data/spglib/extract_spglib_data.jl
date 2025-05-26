using Spglib
using SpaceGroups


function get_elements(hall_number::Integer)::Vector{SpaceGroupElement{3,Int}}
    r, t=get_symmetry_from_database(hall_number)
    n=length(r)
    elements=Vector{SpaceGroupElement{3,Int}}(undef, n)
    for i=1:n
        elements[i]=SpaceGroupElement{3,Int}(r[i], rationalize.(t[i], tol=1e-2))
    end
    elements
end

function expand(G::SpaceGroupQuotient{N, T}, g::SpaceGroupElement{N, T}) where {N, T<:Integer}
    SpaceGroupQuotient{N, T}(union(G.e, Set([g])))
end

function greedy_generators(elements::Vector{SpaceGroupElement{3,Int}})::Vector{SpaceGroupElement{3,Int}}
    generators=Vector{SpaceGroupElement{3,Int}}()
    G=SpaceGroupQuotient{3,Int}() # Start with the trivial group
    candidates=Set(elements)
    while !isempty(candidates)
        # Find the candidate which expands the group the most
        best_candidate=elements[1]
        best_size=0
        for candidate in candidates
            G1=expand(G, candidate)
            size=length(G1.e)
            if size>best_size
                best_size=size
                best_candidate=candidate
            end
        end
        # Add the best candidate to the generators
        push!(generators, best_candidate)
        # expand the group
        G=expand(G, best_candidate)
        # Remove the expanded group elements from the candidates
        setdiff!(candidates, G.e)
    end
    if length(G)>length(elements)
        error("Generated group is larger than the original group.")
    end
    generators
end

function sge_to_tuple_string(g::SpaceGroupElement{3,Int})::String
    t=join(string.(g.b), ", ")
    "($(g.a), [$t])"
end