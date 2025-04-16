# Two-dimensional space groups (wallpaper groups) 
using StaticArrays

#  Rotations
r2=[-1 0; 0 -1] # Two-fold rotation, or central symmetry
r3=[-1 -1; 1 0] # Three-fold rotation
r4=[0 -1; 1 0] # Four-fold rotation
r6=[0 -1; 1 1] # Six-fold rotation

s2=@SGE(r2) # Two-fold rotation
s3=@SGE(r3) # Three-fold rotation
s4=@SGE(r4) # Four-fold rotation
s6=@SGE(r6) # Six-fold rotation

# Mirror operations
m1=[-1 0; 0 1] 
m2=[1 0; 0 -1] 
m3=[0 1; 1 0]
m4=[0 -1; -1 0]

sm1=@SGE(m1)
sm2=@SGE(m2)
sm3=@SGE(m3)
sm4=@SGE(m4)

# Translations
t1=[0, 1//2] 
t2=[1//2, 1//2] 

# Glide reflections
g1=@SGE(m1, t1) 
g2=@SGE(m1, t2) 
g3=@SGE(m3, t2) 

# Space groups
p1=SpaceGroupQuotient{2,Int}() # P1
p2=SpaceGroupQuotient([s2]) # P2
p3=SpaceGroupQuotient([s3]) # P3
p4=SpaceGroupQuotient([s4]) # P4
p6=SpaceGroupQuotient([s6]) # P6

pm=SpaceGroupQuotient([sm1]) # PM
pg=SpaceGroupQuotient([g1]) # PG
cm=SpaceGroupQuotient([sm3]) # CM
pmm=SpaceGroupQuotient([sm1, sm2]) # Pmm
pmg=SpaceGroupQuotient([sm2, g1]) # Pmg

pgg=SpaceGroupQuotient([s2, g2]) # Pgg
cmm=SpaceGroupQuotient([sm3, sm4]) # Cmm
p4m=SpaceGroupQuotient([s4, sm1]) # P4m
p4g=SpaceGroupQuotient([s4, g3]) # P4g
p3m1=SpaceGroupQuotient([s3, sm3]) # P3m1

p31m=SpaceGroupQuotient([s3, sm4]) # P31m
p6m=SpaceGroupQuotient([s6, sm3]) # P6m