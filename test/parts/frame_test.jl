@testset "create frame" begin
    @test_nowarn EulerParameterFrame([0.0,0,0], [1.0, 0, 0, 0])
end
