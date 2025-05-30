using Spglib, JSON
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

function get_group_data(hall_number::Integer)::Dict{String, Any}
    elements=get_elements(hall_number)
    generators=greedy_generators(elements)
    sgt=get_spacegroup_type(hall_number)
    fields=fieldnames(typeof(sgt))
    group_data=Dict{String, Any}()
    for field in fields
        group_data[string(field)]=getfield(sgt, field)
    end
    group_data["generators"]=[sge_to_tuple_string(g) for g in generators]
    group_data["order"]=length(elements)
    group_data
end

function save_groups3D()
    data=Vector{Any}()
    for hall_number in 1:530
        gdata=get_group_data(hall_number)
        push!(data, gdata)
    end
    # Save the data to a JSON file
    script_dir = @__DIR__
    json_file = joinpath(script_dir, "groups3D.json")
    open(json_file, "w") do io
        JSON.print(io, data, 4)
    end
end