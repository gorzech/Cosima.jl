using Cosima
import Cosima: approximate_jacobian, q_mul, rot, q_axis, set_body_coordinates!
using Test
using LinearAlgebra
using StaticArrays

@testset "Cosima.jl" begin
    include("helpers/matrix_test.jl")
    include("parts/frame_test.jl")
    include("parts/rigid_test.jl")
    include("interactions/forces_test.jl")
end
