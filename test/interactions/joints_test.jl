
function approx_constr_derivatives(sys, q, h)
    """Helper test function to approximate contraint derivatives"""
    function normalize_Euler_params(sys)
        for b in sys.bodies.r_bodies
            normalize!(@view b.q[4:7])
        end
        nothing
    end

    function constraint_function(qc_)
        set_body_coordinates!(sys.bodies, qc_)
        # Just in case
        normalize_Euler_params(sys)
        c, _, _, _ = constraints(sys)
        return c
    end

    c_q = approximate_jacobian(constraint_function, q)

    # Do this when they will depend on time
    # # c_p = c_q * h + c_t
    # def constr_t_fun(t):
    #     c, _, _, _ = sys.constraints(t, q, h)
    #     return c

    # c_t = approximate_jacobian(constr_t_fun, t_base)

    function gamma_function(qc_)
        set_body_coordinates!(sys.bodies, qc_)
        normalize_Euler_params(sys)
        _, b_cq, _, _ = constraints(sys)
        g = b_cq * h  # Cq * qp
        return g
    end

    c_p = gamma_function(q) #+ c_t

    cq_qp_dq = approximate_jacobian(gamma_function, q)

    function transf_q_h!(body, q, A_q, A_h)
        A_h[:, body.hi[1:3]] .= A_q[:, body.qi[1:3]]
        LiT = Cosima.gep(q[body.qi[4:7]])'
        A_h[:, body.hi[4:6]] .= 0.5 .* A_q[:, body.qi[4:7]] * LiT
        A_h[:, body.hi[7:end]] .= A_q[:, body.qi[8:end]]
    end

    c_h = zeros(size(c_q, 1), length(h))
    ch_h_dh = zeros(size(c_h))
    # transform _dq terms to correspond _dh
    for b in sys.bodies.r_bodies
        transf_q_h!(b, q, c_q, c_h)
        transf_q_h!(b, q, cq_qp_dq, ch_h_dh)
    end

    # # Approximate c_tt using second order formula (after wiki)
    # step = sqrt(eps(1.0)) * 1000
    # c_tt = (
    #     constr_t_fun(t_base + step)
    #     - 2 * constr_t_fun(t_base)
    #     + constr_t_fun(t_base - step)
    # ) / step ** 2

    g = -ch_h_dh * h #- c_tt

    return (c_h, c_p, g)
end

@testset "joint_point_2_rigid_initial" begin
    b = Bodies()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(b, 1.0, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2; p2])

    s0 = [9.0, 2, -1]
    joint = @test_nowarn JointPoint(rb1, rb2, s0)

    sys = Mbs(b, [joint], Force[])

    @test sys.nconstr == 3
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 3

    @test norm(c) < 2e-14
end


@testset "joint_point_2_rigid_rotate_translate" begin
    b = Bodies()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(b, 1.0, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2; p2])

    s0 = [9.0, 2, -1]
    joint = JointPoint(rb1, rb2, s0)

    sys = Mbs(b, [joint], Force[])

    rt = [17.0, -11.2, pi / 4]
    pt = q_axis(1.1324, z)
    At = rot(pt)

    q = [At * (r1 + rt); q_mul(pt, p1); At * (r2 + rt); q_mul(pt, p2)]
    set_body_coordinates!(b, q)

    c, _, _, _ = constraints(sys)
    @test norm(c) < 2e-14

    # Now rotate and translate first body
    q = [At * (r1 + rt); q_mul(pt, p1); r2; p2]
    set_body_coordinates!(b, q)

    # Location of the joint point for the first body:
    st = At * (s0 + rt)

    c_expected = st - s0
    c, _, _, _ = constraints(sys)
    @test c ≈ c_expected rtol = 1e-14 atol = 1e-14
end

@testset "joint_point_4_bodies_approximate" begin
    Random.seed!(14)
    b = Bodies()
    r1_r = [0.0, 1, 0]
    p1_r = normalize(rand(4))
    rb1 = RBody!(b, 1, ones(3), [r1_r; p1_r])
    r2_r = [0.1, 0.2, 0.7]
    p2_r = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2_r; p2_r])

    joints = [
        JointPoint(
            rb1, rb2, rand(3)
        ),  # For rigid bodies this can be anywhere!
    ]

    sys = Mbs(b, joints, Force[])

    @test sys.nconstr == 3
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 3
    
    @test norm(c) < 1e-14
    
    y0 = initial_position(sys)
    h = rand(sys.nh) * 0.7623
    q = y0[1:sys.nq] + rand(sys.nq) * 1e-2
    # normalize quaternions
    for b in sys.bodies.r_bodies
        normalize!(@view q[b.qi[4:7]])
    end
    set_body_coordinates!(b, q, h)

    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) / norm(c_q) < 2e-5
    @test norm(c_p - expected_c_p) < 2e-14
    @test norm(g - expected_g) / norm(g) < 1e-5
end

@testset "joint_simple_rigid_initial" begin
    b = Bodies()
    r1 = [0.0, 0, 0]
    p1 = [1.0, 0, 0, 0]
    rb1 = RBody!(b, 2, [0.1, 1e-3, 1e-3], [r1; p1])

    joint = @test_nowarn JointSimple(rb1)
    @test length(joint.qi) == 6
    @test length(joint.hi) == 6
    @test size(joint.G0_2) == (3, 3)

    sys = Mbs(b, [joint], Force[])

    @test FlexMbd.nc(sys.joints[1]) == 6
    @test sys.nconstr == 6
    @test length(sys.joints) == 1
    
    c, c_q, c_p, gam = constraints(sys)

    @test length(c) == 6
    @test norm(c) < 1e-16
    @test norm(c_p) < 1e-16
    @test norm(gam) < 1e-16
    @test c_q ≈ diagm([1.0, 1, 1, 0.5, 0.5, 0.5]) rtol = 1e-16 atol = 1e-16
end

@testset "joint_simple_two_bodies" begin
    b = Bodies()
    p2 = normalize(rand(4))

    r1 = [0.1, 0.3, 0.2]
    rb1 = RBody!(b, 2, [0.1, 1e-3, 1e-3], [r1; p2])

    fem2 = BeamFem(Beam3d(0.5), 2)
    rm2 = ReducedModelSimple(fem2)
    r2 = [0.0, 0, 0]
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2; -p2])

    joints = [JointSimple(rb1), JointSimple(fb2)]

    sys = Mbs(b, joints, Force[])

    @test sys.nconstr == 12
    @test length(sys.joints) == 2
    @test FlexMbd.nc(sys.joints[1]) == 6
    @test FlexMbd.nc(sys.joints[2]) == 6

    c, _, _, _ = constraints(sys)

    @test length(c) == 12
    @test norm(c) < 1e-16

    rb1.q[2] = 0.0
    c, _, _, _ = constraints(sys)
    @test abs(norm(c) - 0.3) < 1e-16
end

@testset "joint_simple_two_bodies_approx" begin
    b = Bodies()
    p2 = normalize(rand(4))

    r1 = [0.1, 0.3, 0.2]
    rb1 = RBody!(b, 2, [0.1, 1e-3, 1e-3], [r1; p2])

    fem2 = BeamFem(Beam3d(0.5), 2)
    rm2 = ReducedModelSimple(fem2)
    r2 = [0, 0, 0]
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2; -p2])

    joints = [JointSimple(rb1), JointSimple(fb2)]

    sys = Mbs(b, joints, Force[])

    @test sys.nconstr == 12
    @test length(sys.joints) == 2

    q0 = initial_position(sys)
    q = q0[1:sys.nq] .+ 1e-6
    h = rand(sys.nh) .* 0.7623
    # normalize quaternions
    normalize!(@view q[4:7])
    normalize!(@view q[11:14])
    set_body_coordinates!(b, q, h)
    
    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) < 1e-5
    @test norm(c_p - expected_c_p) < 1e-14
    @test norm(g - expected_g) < 1e-14
end

@testset "joint_simple_approx_not_all_rigid" begin
    b = Bodies()
    p2 = normalize(rand(4))

    r1 = [0.1, 0.3, 0.2]
    rb1 = RBody!(b, 2, [0.1, 1e-3, 1e-3], [r1; p2])

    fem2 = BeamFem(Beam3d(0.5), 2)
    rm2 = ReducedModelSimple(fem2)
    r2 = [0.0, 0, 0]
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2; -p2])

    joints = [JointSimple(rb1, [1, 2], false), JointSimple(fb2, Int[])]

    sys = Mbs(b, joints, Force[])

    @test sys.nconstr == 5
    @test length(sys.joints) == 2
    @test FlexMbd.nc(sys.joints[1]) == 2
    @test FlexMbd.nc(sys.joints[2]) == 3

    q0 = initial_position(sys)
    q = q0[1:sys.nq] .+ 1e-6
    h = rand(sys.nh) .* 0.7623
    # normalize quaternions
    normalize!(@view q[4:7])
    normalize!(@view q[11:14])
    set_body_coordinates!(b, q, h)
    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) < 8e-6
    @test norm(c_p - expected_c_p) < 1e-14
    @test norm(g - expected_g) < 1e-14
end

@testset "joint_perped1_2_rigid_initial" begin
    b = Bodies()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(b, 1, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2; p2])

    v_1 = SA[1.0, 0, 0]
    v_2 = SA[0.0, 1, 0]
    joint = @test_nowarn JointPerpend1(rb1, rb2, v_1, v_2)

    sys = Mbs(b, [joint], Force[])

    @test sys.nconstr == 1
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 1

    @test abs(c[1]) < 1e-14
end

@testset "joint_perped1_2_rigid_rotate_translate" begin
    b = Bodies()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(b, 1, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2; p2])

    v_1 = [1.0, 0, 0]
    v_2 = [0.0, 1, 0]
    joint = JointPerpend1(rb1, rb2, v_1, v_2)
    sys = Mbs(b, [joint], Force[])

    rt = [17, -11.2, pi / 4]
    pt = q_axis(1.1324, z)
    At = rot(pt)

    q = [At * (r1 + rt); q_mul(pt, p1); At * (r2 + rt); q_mul(pt, p2)]
    set_body_coordinates!(b, q)

    c, _, _, _ = constraints(sys)
    @test abs(c[1]) < 1e-14

    # Now rotate and translate first body
    q = [At * (r1 + rt); q_mul(pt, p1); r2; p2]
    set_body_coordinates!(b, q)

    # Location of the joint point for the first body:
    v_1_n = At * v_1

    c_expected = [v_1_n' * v_2]
    c, _, _, _ = constraints(sys)
    @test c ≈ c_expected rtol = 1e-14 atol = 1e-14
end

@testset "joint_perped1_2_flexible_initial" begin
    b = Bodies()
    s0 = [9.0, 2, -1]
    s1_local = [4.5, 0, 0]
    p1 = normalize(rand(4))
    r1 = s0 - rot(p1) * s1_local
    fem1 = BeamFem(Beam3d(6.0 / 4), 4)
    rm1 = ReducedModelSimple(fem1, at_center_of_mass=false)
    fb1 = FBodyDirectIntegration!(b, fem1, rm1, [r1; p1])

    p2 = normalize(rand(4))
    fem2 = BeamFem(Beam3d(0.2 / 3), 3)
    rm2 = ReducedModelSimple(fem2, at_center_of_mass=false)
    s2_local = [rm2.qn0[6], 0, 0]
    r2 = s0 - rot(p2) * s2_local
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2; p2])

    v_1 = [1.0, 0, 0]
    v_2 = [0.0, 1, 0]
    # For flexible bodies nonzero locations are needed
    joint = JointPerpend1(fb1, fb2, v_1, v_2, s0)

    sys = Mbs(b, [joint], Force[])

    @test sys.nconstr == 1
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 1

    @test abs(c[1]) < 1e-14
end

@testset "joint_perped1_2_flexible_rotate_translate" begin
    b = Bodies()
    s0 = [9.0, 2, -1]
    s1_local = [4.5, 0, 0]
    p1 = normalize(rand(4))
    r1 = s0 - rot(p1) * s1_local
    fem1 = BeamFem(Beam3d(6.0 / 4), 4)
    rm1 = ReducedModelSimple(fem1, at_center_of_mass=false)
    fb1 = FBodyDirectIntegration!(b, fem1, rm1, [r1; p1])

    p2 = normalize(rand(4))
    fem2 = BeamFem(Beam3d(0.2 / 3), 3)
    rm2 = ReducedModelSimple(fem2, at_center_of_mass=false)
    s2_local = [rm2.qn0[6], 0, 0]
    r2 = s0 - rot(p2) * s2_local
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2; p2])

    v_1 = SA[1.0, 0, 0]
    v_2 = SA[0.0, 1, 0]
    # For flexible bodies nonzero locations are needed
    joint = JointPerpend1(fb1, fb2, v_1, v_2, s0)

    sys = Mbs(b, [joint], Force[])

    rt = [17, -11.2, pi / 4]
    pt = q_axis(1.1324, z)
    At = rot(pt)

    q = [
        At * (r1 + rt);
        q_mul(pt, p1);
        zeros(fb1.nflex);
        At * (r2 + rt);
        q_mul(pt, p2);
        zeros(fb2.nflex)
    ]
    set_body_coordinates!(b, q)

    c, _, _, _ = constraints(sys)
    @test abs(c[1]) < 1e-14

    # Now rotate and translate second body
    q = [
        r1;
        p1;
        zeros(fb1.nflex);
        At * (r2 + rt);
        q_mul(pt, p2);
        zeros(fb2.nflex)
    ]
    set_body_coordinates!(b, q)

    # Location of the joint point for the second body:
    v_2_n = At * v_2

    c_expected = v_1' * v_2_n
    c, _, _, _ = constraints(sys)
    @test c[1] ≈ c_expected rtol = 1e-14 atol = 1e-14

    # Add some deformation to the flexible body (first)
    q = initial_position(sys)
    q[23:25] .= 0.01
    set_body_coordinates!(b, q)

    At = skew(q[23:25]) + I
    # Unit vector at the same body
    A1 = rot(p1)
    v_1_n = A1 * At * A1' * v_1

    c_expected = [v_1_n' * v_2]
    c, _, _, _ = constraints(sys)
    @test c ≈ c_expected rtol = 1e-14 atol = 1e-14
end

@testset "joint_perped1_4_bodies_approximate" begin
    Random.seed!(856)
    b = Bodies()
    r1_r = [0.0, 1, 0]
    p1_r = normalize(rand(4))
    rb1 = RBody!(b, 1, ones(3), [r1_r; p1_r])
    r2_r = [0.1, 0.2, 0.7]
    p2_r = normalize(rand(4))
    rb2 = RBody!(b, 1, ones(3), [r2_r; p2_r])

    s0 = [-9.0, 12, 1 / 3]
    s1_local = [4.5, 0, 0]
    p1_f = normalize(rand(4))
    r1_f = s0 - rot(p1_f) * s1_local
    fem1 = BeamFem(Beam3d(6.0 / 4), 4)
    rm1 = ReducedModelSimple(fem1, at_center_of_mass=false)
    fb1 = FBodyDirectIntegration!(b, fem1, rm1, [r1_f; p1_f])

    p2_f = normalize(rand(4))
    fem2 = BeamFem(Beam3d(0.2 / 3), 3)
    rm2 = ReducedModelSimple(fem2, at_center_of_mass=false)
    s2_local = [rm2.qn0[6], 0, 0]
    r2_f = s0 - rot(p2_f) * s2_local
    fb2 = FBodyDirectIntegration!(b, fem2, rm2, [r2_f; p2_f])

    v_1 = [1.0, 0, 0]
    v_2 = [0, 0.3, 0.4]
    # For flexible bodies locations are needed

    joints = [
        JointPerpend1(rb1, rb2, v_1, v_2),  # For rigid bodies this can be anywhere!
        JointPerpend1(fb2, rb2, v_1, [0, -0.1, -0.07], s0),
        JointPerpend1(fb1, fb2, [0.0, 7, 3], [0.0, -3, 7], s0),
    ]

    sys = Mbs(b, joints, Force[])

    @test sys.nconstr == 3
    @test length(sys.joints) == 3

    c, _, _, _ = constraints(sys)
    @test length(c) == 3
    @test norm(c) < 1e-14

    y0 = initial_position(sys)
    h = rand(sys.nh) .* 0.7623
    q = y0[1:sys.nq] .+ rand(sys.nq) * 1e-2
    # normalize quaternions
    for b in sys.bodies.r_bodies
        normalize!(@view q[b.qi[4:7]])
    end
    set_body_coordinates!(b, q, h)

    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) / norm(c_q) < 8e-6
    @test norm(c_p - expected_c_p) < 2e-14
    @test norm(g - expected_g) < 3e-5
end