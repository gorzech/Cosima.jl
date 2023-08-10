module Cosima
using StaticArrays
# Write your package code here.

export EulerParameterFrame

include("parts/frame.jl")

include("parts/rigid.jl")

end
