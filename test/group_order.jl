# Test the group order of various space groups

include("util/icosahedral.jl")

@testset "Icosahedral groups" begin
        @test length(PI)==60
        @test length(PIh)==120
        @test length(PI_n)==60
        @test length(PIh_n)==120
end

