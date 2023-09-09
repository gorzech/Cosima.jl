module Cosima
using StaticArrays
using LinearAlgebra
# Write your package code here.

export EulerParameterFrame

include("helpers/matrix.jl")
export x, y, z

include("parts/node.jl")

include("parts/rigid.jl")
export RBody!, mass, Bodies

include("interactions/forces.jl")
export ForceTorque, Force

include("interactions/joints.jl")
export Joint, JointSimple, JointPoint, JointPerpend1

include("globals/system.jl")
export Mbs, force, OdeMbs, initial_position, ode!, constraints

end
