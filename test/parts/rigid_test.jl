@testset "mass_single_rigid_body" begin
    m = 12.1
    R3 = rand(3, 3)
    Ic = 0.5 * (R3 + R3') # make symmetric
    sys = Mbs()
    b = RBody!(sys, m, Ic)

    mass_test = mass(b)

    mass_expected = diagm([m, m, m, 0, 0, 0])
    mass_expected[4:end, 4:end] = Ic
    @test mass_test == mass_expected
end

@testset "grav_single_rigid_body" begin
    g = [0.1, -9.81, 6]

    m = 3.2
    Ic = [1.2, 0.13, 78e-3]
    sys = Mbs()
    b = RBody!(sys, m, Ic)

    expected_g = [m * g; zeros(3)]

    test_g = Cosima.grav(b, g)

    @test test_g ≈ expected_g
end

@testset "quadratic_inertia_single_rigid_body" begin
    Ic = diagm([1.2, 0.13, 78e-3])
    sys = Mbs()
    b = RBody!(sys, 3.2, Ic)

    om = rand(3)
    sys.bodies.r_bodies[1].node.h[4:end] .= om

    expected_b = [zeros(3); Cosima.skew(om) * Ic * om]

    test_b = Cosima.quadratic_inertia(b)

    @test expected_b ≈ test_b rtol = 1e-12 atol = 1e-14
end