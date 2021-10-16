module SpaceGroups

import StaticArrays: SMatrix,SVector

export SpaceGroup, SpaceGroupElement

include("finitegroup.jl")
include("spacegroup.jl")

end
