@testset "mass_single_body" begin
    m = 12.1
    R3 = rand(3, 3)
    Ic = 0.5 * (R3 + R3') # make symmetric
    bodies = Bodies()
    b = RBody!(bodies, m, Ic)
    sys = Mbs(bodies, Joint[], Force[])

    mass_test = mass(sys)

    mass_expected = diagm([m, m, m, 0, 0, 0])
    mass_expected[4:end, 4:end] = Ic
    @test mass_test == mass_expected
end

@testset "mass_three_rigid_bodies" begin
    bodies = Bodies()
    RBody!(bodies, 1.0, ones(3) .* 1.2),
    RBody!(bodies, 2.0, ones(3) .* 3.7),
    RBody!(bodies, 3.0, ones(3) .* 7.9)

    sys = Mbs(bodies, Joint[], Force[])

    mass_expected = diagm(repeat([1.0 1.2 2 3.7 3 7.9], 3, 1)[:])
    mass_test = mass(sys)
    @test mass_test == mass_expected
end

@testset "grav_two_bodies" begin
    g = SA[0.1, -9.81, 6]

    m1 = 3.2
    m2 = 2
    Ic1 = rand(3, 3)
    Ic2 = [1.2, 0.13, 78e-3]
    bodies = Bodies()
    RBody!(bodies, m1, Ic1), RBody!(bodies, m2, Ic2)
    sys = Mbs(bodies, Joint[], Force[], gv=g, use_quadratic=false)

    expected_g = [m1 * g; zeros(3); m2 * g; zeros(3)]

    test_g = force(0.0, sys)

    @test test_g ≈ expected_g
end

@testset "quadratic_inertia" begin
    Ic1 = rand(3, 3)
    Ic2 = diagm([1.2, 0.13, 78e-3])
    bodies = Bodies()
    RBody!(bodies, 3.2, Ic1), RBody!(bodies, 2, Ic2)

    om1 = rand(3)
    om2 = [0, 1e-2, 0.1234112]
    bodies.r_bodies[1].node.h[4:end] .= om1
    bodies.r_bodies[2].node.h[4:end] .= om2

    sys = Mbs(bodies, Joint[], Force[])

    expected_b = [zeros(3); skew(om1) * Ic1 * om1; zeros(3); skew(om2) * Ic2 * om2]

    test_b = -force(0.0, sys)

    @test expected_b ≈ test_b rtol = 1e-12 atol = 1e-14
end

@testset "ode_free_fall_rigid_single" begin
    bodies = Bodies()
    RBody!(bodies, 1.17893, ones(3), SA[1.0, 0, 2, 1, 0, 0, 0])

    g = 9.81
    grav = [0, 0, -g]
    sys = Mbs(bodies, Joint[], Force[], gv=grav)

    y0 = initial_position(sys)
    prob = ODEProblem(ode!, y0, (0.0, 1.0), OdeMbs(sys))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success
    tend = sol.t[end]
    test_y = sol[end]

    expected_y = copy(y0)
    expected_y[3] -= 0.5 * g * tend^2
    expected_y[10] -= g * tend

    @test expected_y ≈ test_y rtol = 1e-13 atol = 1e-14
end

@testset "ode_rotation_const_angular_velocity" begin
    b = Bodies()
    RBody!(b, 1.17893, ones(3), [1.0, 0, 2, 1, 0, 0, 0])
    sys = Mbs(b, Joint[], Force[])

    y0 = initial_position(sys)

    # set angular velocity
    for ii = 11:13
        y0[ii] = 2 * pi
        expected_y = copy(y0)
        prob = ODEProblem(ode!, y0, (0.0, 1.0), OdeMbs(sys))
        sol = solve(prob, saveat=0.5, reltol=1e-6, abstol=1e-9)

        # res = solve_ivp(sys, (0, 1), y0, t_eval=(0.5, 1.0), max_step=0.07)
        # assert res.success, "Integration fails"
        @test sol.retcode == ReturnCode.Success
        y0[ii] = 0

        # tend = T(end);
        cmp_idx = [1; 2; 3; 8:length(y0)]

        t_idx = 2  # half of the motion - should rotate by 180 deg
        @test abs(sol.t[t_idx] - 0.5) < 1e-15
        @test sol[cmp_idx, t_idx] ≈ expected_y[cmp_idx] rtol = 1e-14 atol = 1e-14

        A_expected = rot_axis(pi, Cosima.axis(ii - 10))  # Rotation by 180 deg
        p = normalize(sol[4:7, t_idx])
        A_idx = rot(p)
        @test A_idx ≈ A_expected rtol = 1e-7 atol = 1e-8#"Not equal at t=0.5: {np.linalg.norm(A_idx - A_expected)}"

        t_idx = 3  # end
        @test abs(sol.t[t_idx] - 1.0) < 1e-15
        @test sol[cmp_idx, t_idx] ≈ expected_y[cmp_idx] rtol = 1e-14 atol = 1e-14
        # As EP are not unique - check if their rotational matrices are the same
        # A_expected = np.eye(3)  # body should perform full rotation
        p = normalize(sol[4:7, t_idx])
        A_idx = rot(p)
        @test A_idx ≈ I rtol = 1e-6 atol = 1e-7 #"Not equal at t=1: {np.linalg.norm(A_idx - A_expected)}"
    end
end

@testset "ode_two_bodies_rotation_translation" begin
    Ic2 = diagm([7.12, 0.23, 0.87])
    p2 = normalize(rand(4))
    om2_0 = [0.1, 0, 0.3]

    b = Bodies()
    RBody!(b, 1.17893, ones(3), [1.0, 0, 2, 1, 0, 0, 0])
    RBody!(b, 289.9872, Ic2, [0.0; 0; 0; p2], [0; 10.0; 0; om2_0])

    g = 9.81
    gv = [0, 0, -g]
    sys = Mbs(b, Joint[], Force[], gv=gv)
    y0 = initial_position(sys)
    prob = ODEProblem(ode!, y0, (0.0, 1.0), OdeMbs(sys))
    sol = solve(prob, saveat=1.0, reltol=1e-6, abstol=1e-9)
    @test sol.retcode == ReturnCode.Success

    t_end = sol.t[2]
    y_end = sol[:, 2]

    y_expected = copy(y0)
    y_expected[3] -= 0.5 * g * t_end^2
    y_expected[9] += 10 * t_end  # due to initial velocity
    y_expected[10] -= 0.5 * g * t_end^2
    y_expected[17] -= g * t_end
    y_expected[23] -= g * t_end
    cmp_idx = [1:10; 15:23]
    @test y_end[cmp_idx] ≈ y_expected[cmp_idx] rtol = 1e-13 atol = 1e-14

    # Check if EP are normalized and angular kinetic energy preserved for body 1
    p2 = y_end[11:14]
    @test abs(p2' * p2 - 1) < 1e-9
    om2 = y_end[24:26]
    Eak_end = 0.5 * om2' * Ic2 * om2
    Eak_0 = 0.5 * om2_0' * Ic2 * om2_0
    @test abs(Eak_0 - Eak_end) < 1e-9 #"Difference in angular kinetic energy too large: $(abs(Eak_0 - Eak_end))"
end