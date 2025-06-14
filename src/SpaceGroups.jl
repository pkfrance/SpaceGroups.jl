module SpaceGroups

import StaticArrays: SMatrix, SVector
import LinearAlgebra: ⋅, I, rank
import Base: *, ∘

export SpaceGroupQuotient, SpaceGroupElement, @SGE, ∘
export WyckoffPosition, stabilizer_quotient, is_valid_wyckoff, normalize
export AffinePhase, ComplexOrbit, RealOrbit, ExtinctOrbit, FormalOrbit, PhysicalOrbit
export make_orbit


include("finitegroup.jl")
include("spacegroup.jl")
include("wyckoff.jl")
include("orbit.jl")

end
