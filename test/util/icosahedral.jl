# Define some space groups for icosahedral quasicrystals

using StaticArrays


############################## Groups of type "P" (primitive) ##############################
# Three-fold rotation
r3 = 
[0 0 1 0 0 0
 1 0 0 0 0 0
 0 1 0 0 0 0
 0 0 0 0 0 1
 0 0 0 1 0 0
 0 0 0 0 1 0]

# Five fold rotation
r5 = 
[1 0 0 0 0 0
 0 0 0 0 1 0
 0 1 0 0 0 0
 0 0 1 0 0 0
 0 0 0 0 0 -1
 0 0 0 -1 0 0]

# Central symmetry
c = 
[-1 0 0 0 0 0;
 0 -1 0 0 0 0;
 0 0 -1 0 0 0;
 0 0 0 -1 0 0;
 0 0 0 0 -1 0;
 0 0 0 0 0 -1]

# Translations for non-symmorphic icosahedral groups
t1 = [1 // 5, 0 // 1, -1 // 5, 0 // 1, 0 // 1, -1 // 5]
t2 = [0 // 1, 0 // 1, 0 // 1, 0 // 1, 1 // 2, -1 // 2]

# Rotations
s3 = @SGE(r3)
s5 = @SGE(r5)

# Central symmetry
sc = @SGE(c)

# Non-symmorphic operations
s5_1 = @SGE(r5, t1)
s5_2 = @SGE(r5, t2)

const PI=SpaceGroupQuotient([s3,s5]) # Symmorphic non-centrosymmetric group
const PIh=SpaceGroupQuotient([s3,s5,sc]) # Symmorphic centrosymmetric group
const PI_n=SpaceGroupQuotient([s3,s5_1]) # Non-symmorphic non-centrosymmetric group
const PIh_n=SpaceGroupQuotient([s3,s5_2,sc]) # Non-symmorphic centrosymmetric group