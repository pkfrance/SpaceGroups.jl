module SpaceGroups

import StaticArrays: SMatrix, SVector
import LinearAlgebra: ⋅

export SpaceGroupQuotient, SpaceGroupElement
export WyckoffPosition
export AffinePhase, ComplexOrbit, RealOrbit, ExtinctOrbit, FormalOrbit, PhysicalOrbit
export make_orbit


include("finitegroup.jl")
include("spacegroup.jl")
include("wyckoff.jl")
include("orbit.jl")

end
