module Cosima
using StaticArrays
using LinearAlgebra
# Write your package code here.

export EulerParameterFrame

include("parts/frame.jl")

include("parts/rigid.jl")
export RBody!, mass, Bodies

end
