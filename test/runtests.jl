using Cosima
import Cosima: approximate_jacobian, q_mul, rot, q_axis, set_body_coordinates!, skew, rot_axis
using Test
using LinearAlgebra
using StaticArrays
using Random
using DifferentialEquations

@testset "Cosima.jl" begin
    include("helpers/matrix_test.jl")
    include("parts/rigid_test.jl")
    include("interactions/forces_test.jl")
    include("interactions/joints_test.jl")
    include("globals/system_test.jl")
end
