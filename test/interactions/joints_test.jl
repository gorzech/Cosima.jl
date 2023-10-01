
function approx_constr_derivatives(sys, q, h)
    """Helper test function to approximate contraint derivatives"""
    function normalize_Euler_params(sys)
        for b in sys.bodies.r_bodies
            normalize!(@view b.node.q[4:7])
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
        A_h[:, body.node.hi[1:3]] .= A_q[:, body.node.qi[1:3]]
        LiT = Cosima.gep(q[body.node.qi[4:7]])'
        A_h[:, body.node.hi[4:6]] .= 0.5 .* A_q[:, body.node.qi[4:7]] * LiT
        A_h[:, body.node.hi[7:end]] .= A_q[:, body.node.qi[8:end]]
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
    sys = Mbs()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(sys, 1.0, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2; p2])

    s0 = [9.0, 2, -1]
    joint = @test_nowarn JointPoint!(sys, rb1, rb2, s0)

    @test sys.nconstr == 3
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 3

    @test norm(c) < 2e-14
end


@testset "joint_point_2_rigid_rotate_translate" begin
    sys = Mbs()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(sys, 1.0, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2; p2])

    s0 = [9.0, 2, -1]
    joint = JointPoint!(sys, rb1, rb2, s0)

    rt = [17.0, -11.2, pi / 4]
    pt = q_axis(1.1324, z)
    At = rot(pt)

    q = [At * (r1 + rt); q_mul(pt, p1); At * (r2 + rt); q_mul(pt, p2)]
    set_body_coordinates!(sys.bodies, q)

    c, _, _, _ = constraints(sys)
    @test norm(c) < 2e-14

    # Now rotate and translate first body
    q = [At * (r1 + rt); q_mul(pt, p1); r2; p2]
    set_body_coordinates!(sys.bodies, q)

    # Location of the joint point for the first body:
    st = At * (s0 + rt)

    c_expected = st - s0
    c, _, _, _ = constraints(sys)
    @test c ≈ c_expected rtol = 1e-14 atol = 1e-14
end

@testset "joint_point_4_bodies_approximate" begin
    Random.seed!(14)
    sys = Mbs()
    r1_r = [0.0, 1, 0]
    p1_r = normalize(rand(4))
    rb1 = RBody!(sys, 1, ones(3), [r1_r; p1_r])
    r2_r = [0.1, 0.2, 0.7]
    p2_r = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2_r; p2_r])

    JointPoint!(
        sys, rb1, rb2, rand(3)
    )

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
        normalize!(@view q[b.node.qi[4:7]])
    end
    set_body_coordinates!(sys.bodies, q, h)

    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) / norm(c_q) < 2e-5
    @test norm(c_p - expected_c_p) < 2e-14
    @test norm(g - expected_g) / norm(g) < 1e-5
end

@testset "joint_simple_rigid_initial" begin
    sys = Mbs()
    r1 = [0.0, 0, 0]
    p1 = [1.0, 0, 0, 0]
    rb1 = RBody!(sys, 2, [0.1, 1e-3, 1e-3], [r1; p1])

    joint = @test_nowarn JointSimple!(sys, rb1)
    @test length(joint.qi) == 6
    @test length(joint.hi) == 6
    @test size(joint.G0_2) == (3, 3)

    @test Cosima.nc(sys.joints[1]) == 6
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
    sys = Mbs()
    p2 = normalize(rand(4))

    r1 = [0.1, 0.3, 0.2]
    rb1 = RBody!(sys, 2, [0.1, 1e-3, 1e-3], [r1; p2])

    r2 = [0.0, 0, 0]
    rb2 = RBody!(sys, 2, [0.1, 1e-3, 1e-3], [r2; -p2])

    JointSimple!(sys, rb1)
    JointSimple!(sys, rb2)

    @test sys.nconstr == 12
    @test length(sys.joints) == 2
    @test Cosima.nc(sys.joints[1]) == 6
    @test Cosima.nc(sys.joints[2]) == 6

    c, _, _, _ = constraints(sys)

    @test length(c) == 12
    @test norm(c) < 1e-16

    rb1.node.q[2] = 0.0
    c, _, _, _ = constraints(sys)
    @test abs(norm(c) - 0.3) < 1e-16
end

@testset "joint_simple_two_bodies_approx" begin
    sys = Mbs()
    p2 = normalize(rand(4))

    r1 = [0.1, 0.3, 0.2]
    rb1 = RBody!(sys, 2, [0.1, 1e-3, 1e-3], [r1; p2])

    r2 = [0.0, 0, 0]
    rb2 = RBody!(sys, 2, [0.1, 1e-3, 1e-3], [r2; -p2])

    JointSimple!(sys, rb1)
    JointSimple!(sys, rb2)

    @test sys.nconstr == 12
    @test length(sys.joints) == 2

    q0 = initial_position(sys)
    q = q0[1:sys.nq] .+ 1e-6
    h = rand(sys.nh) .* 0.7623
    # normalize quaternions
    normalize!(@view q[4:7])
    normalize!(@view q[11:14])
    set_body_coordinates!(sys.bodies, q, h)

    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) < 1e-5
    @test norm(c_p - expected_c_p) < 1e-14
    @test norm(g - expected_g) < 1e-14
end

@testset "joint_perped1_2_rigid_initial" begin
    sys = Mbs()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(sys, 1, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2; p2])

    v_1 = SA[1.0, 0, 0]
    v_2 = SA[0.0, 1, 0]
    joint = @test_nowarn JointPerpend1!(sys, rb1, rb2, v_1, v_2)

    @test sys.nconstr == 1
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 1

    @test abs(c[1]) < 1e-14
end

@testset "joint_perped1_2_rigid_rotate_translate" begin
    sys = Mbs()
    r1 = [0.0, 1, 0]
    p1 = normalize(rand(4))
    rb1 = RBody!(sys, 1, ones(3), [r1; p1])
    r2 = [0.1, 0.2, 0.7]
    p2 = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2; p2])

    v_1 = [1.0, 0, 0]
    v_2 = [0.0, 1, 0]
    joint = JointPerpend1!(sys, rb1, rb2, v_1, v_2)

    rt = [17, -11.2, pi / 4]
    pt = q_axis(1.1324, z)
    At = rot(pt)

    q = [At * (r1 + rt); q_mul(pt, p1); At * (r2 + rt); q_mul(pt, p2)]
    set_body_coordinates!(sys.bodies, q)

    c, _, _, _ = constraints(sys)
    @test abs(c[1]) < 1e-14

    # Now rotate and translate first body
    q = [At * (r1 + rt); q_mul(pt, p1); r2; p2]
    set_body_coordinates!(sys.bodies, q)

    # Location of the joint point for the first body:
    v_1_n = At * v_1

    c_expected = [v_1_n' * v_2]
    c, _, _, _ = constraints(sys)
    @test c ≈ c_expected rtol = 1e-14 atol = 1e-14
end

@testset "joint_perped1_4_bodies_approximate" begin
    Random.seed!(856)
    sys = Mbs()
    r1_r = [0.0, 1, 0]
    p1_r = normalize(rand(4))
    rb1 = RBody!(sys, 1, ones(3), [r1_r; p1_r])
    r2_r = [0.1, 0.2, 0.7]
    p2_r = normalize(rand(4))
    rb2 = RBody!(sys, 1, ones(3), [r2_r; p2_r])

    v_1 = [1.0, 0, 0]
    v_2 = [0, 0.3, 0.4]
    # For flexible bodies locations are needed

    JointPerpend1!(sys, rb1, rb2, v_1, v_2),  # For rigid bodies this can be anywhere!
    @test sys.nconstr == 1
    @test length(sys.joints) == 1

    c, _, _, _ = constraints(sys)
    @test length(c) == 1
    @test norm(c) < 1e-14

    y0 = initial_position(sys)
    h = rand(sys.nh) .* 0.7623
    q = y0[1:sys.nq] .+ rand(sys.nq) * 1e-2
    # normalize quaternions
    for b in sys.bodies.r_bodies
        normalize!(@view q[b.node.qi[4:7]])
    end
    set_body_coordinates!(sys.bodies, q, h)

    _, c_q, c_p, g = constraints(sys)
    expected_c_h, expected_c_p, expected_g = approx_constr_derivatives(sys, q, h)

    @test norm(c_q - expected_c_h) / norm(c_q) < 8e-6
    @test norm(c_p - expected_c_p) < 2e-14
    @test norm(g - expected_g) < 3e-5
end