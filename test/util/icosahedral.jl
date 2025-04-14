# Define some space groups for icosahedral quasicrystals

using StaticArrays


############################## Groups of type "P" (primitive) ##############################
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

const PI=SpaceGroupQuotient([s3,s5]) # Symmorphic non-centrosymmetric group
const PIh=SpaceGroupQuotient([s3,s5,sc]) # Symmorphic centrosymmetric group
const PI_n=SpaceGroupQuotient([s3,s5_1]) # Non-symmorphic non-centrosymmetric group
const PIh_n=SpaceGroupQuotient([s3,s5_2,sc]) # Non-symmorphic centrosymmetric group