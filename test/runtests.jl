using Cosima
using Test
using LinearAlgebra
using StaticArrays

@testset "Cosima.jl" begin
    include("parts/frame_test.jl")
    include("parts/rigid_test.jl")
    include("interactions/forces_test.jl")
end
