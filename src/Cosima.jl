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
export Mbs

include("globals/system_equations.jl")
export force, initial_position, constraints

include("globals/system_ode.jl")
export OdeMbs, ode!

end
