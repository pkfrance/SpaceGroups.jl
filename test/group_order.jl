# Test the group order of various space groups

include("util/icosahedral.jl")
include("util/wallpaper.jl")

@testset "Icosahedral groups" begin
        @test length(PI)==60
        @test length(PIh)==120
        @test length(PI_n)==60
        @test length(PIh_n)==120
end

@testset "Wallpaper groups" begin
        @test length(p1)==1
        @test length(p2)==2
        @test length(p3)==3
        @test length(p4)==4
        @test length(p6)==6
        @test length(pm)==2
        @test length(pg)==2
        @test length(cm)==2
        @test length(pmm)==4
        @test length(pmg)==4
        @test length(pgg)==4
        @test length(cmm)==4
        @test length(p4m)==8
        @test length(p4g)==8
        @test length(p3m1)==6
        @test length(p31m)==6
        @test length(p6m)==12
end

