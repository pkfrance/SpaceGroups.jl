using StaticArrays

@testset "Icosahedral symmetry" begin
    
# Three-fold rotation
r3 = SMatrix{6,6}(
   [0 0 1 0 0 0;
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 1 0 0;
    0 0 0 0 1 0])

# Five fold rotation
r5 = SMatrix{6,6}(
   [1 0 0 0 0 0;
    0 0 0 0 1 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 -1;
    0 0 0 -1 0 0])

# Central symmetry
c = SMatrix{6,6}(
   [-1 0 0 0 0 0;
    0 -1 0 0 0 0;
    0 0 -1 0 0 0;
    0 0 0 -1 0 0;
    0 0 0 0 -1 0;
    0 0 0 0 0 -1])

# Translations
# Zero translation
z = zeros(SVector{6,Rational{Int}})

# Translations for non-symmorphic icosahedral groups
t1 = SVector{6,Rational{Int}}([1 // 5, 0 // 1, -1 // 5, 0 // 1, 0 // 1, -1 // 5])
t2 = SVector{6,Rational{Int}}([0 // 1, 0 // 1, 0 // 1, 0 // 1, 1 // 2, -1 // 2])

# Rotations
s3 = SpaceGroupElement(r3, z)
s5 = SpaceGroupElement(r5, z)

# Central symmetry
sc = SpaceGroupElement(c, z)

# Non-symmorphic operations
s5_1 = SpaceGroupElement(r5, t1)
s5_2 = SpaceGroupElement(r5, t2)

    @testset "Symmorphic non-centrosymmetric group" begin
        PI=SpaceGroup([s3,s5])
        @test length(PI)==60
    end

    @testset "Symmorphic centrosymmetric group" begin
        PIh=SpaceGroup([s3,s5,sc])
        @test length(PIh)==120
    end

    @testset "Non-symmorphic non-centrosymmetric group" begin
        PI_n=SpaceGroup([s3,s5_1])
        @test length(PI_n)==60
    end

    @testset "Non-symmorphic centrosymmetric group" begin
        PIh_n=SpaceGroup([s3, s5_2, sc])
        @test length(PIh_n)==120
    end    
end