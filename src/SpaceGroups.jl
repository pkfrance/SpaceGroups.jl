module SpaceGroups

import StaticArrays: SMatrix,SVector

export SpaceGroupQuotient, SpaceGroupElement

include("finitegroup.jl")
include("spacegroup.jl")

end
