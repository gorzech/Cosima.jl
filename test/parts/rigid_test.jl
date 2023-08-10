@testset "mass_single_rigid_body" begin
    m = 12.1
    R3 = rand(3, 3)
    Ic = 0.5 * (R3 + R3') # make symmetric
    bodies = Bodies()
    b = RBody!(bodies, m, Ic)

    mass_test = mass(b)

    mass_expected = diagm([m, m, m, 0, 0, 0])
    mass_expected[4:end, 4:end] = Ic
    @test mass_test == mass_expected
end