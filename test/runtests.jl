using Cosima
using Test
using LinearAlgebra

@testset "Cosima.jl" begin
    include("parts/frame_test.jl")
    include("parts/rigid_test.jl")
end
