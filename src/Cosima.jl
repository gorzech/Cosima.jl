module Cosima
using StaticArrays
using LinearAlgebra
# Write your package code here.

export EulerParameterFrame

include("helpers/matrix.jl")

include("parts/frame.jl")

include("parts/rigid.jl")
export RBody!, mass, Bodies

include("interactions/forces.jl")

include("interactions/joints.jl")

include("globals/system.jl")

end
