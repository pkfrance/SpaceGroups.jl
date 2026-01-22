
using Test
using SpaceGroups


@testset "SpaceGroups.jl" begin
    include("group_order.jl")
    include("smith_normal_form.jl")
end
