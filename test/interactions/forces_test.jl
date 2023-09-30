@testset "torque_single_rigid_body" begin    
    r1 = [0.0, 1, -1]
    p1 = normalize(rand(4))
    sys = Mbs(use_quadratic = false)
    rb1 = RBody!(sys, 1.0, ones(3), [r1; p1])
    
    # torque_point = np.array([0, 1, -1])
    # s0 = r1 + torque_point
    torque(t) = SA[3, 2 + 2 * t, t ^ 2 - 1]
    trq1 = ForceTorque!(sys, rb1, torque)
    
    @test length(sys.forces) == 1
    
    # y0 = initial_position(sys)
    Q = force(0.1, sys)
    @test length(Q) == 6
    
    T = torque(0.1)
    Q_expected = [zeros(3); Cosima.rot(p1)' * T]
    
    @test norm(Q - Q_expected) < 1e-14
end

@testset "torque_support_rigid_body" begin
    m = 1.25
    l = 0.7
    g = 9.81

    # horizontal beam under gravity in y
    # equivalent torque
    torque(_) = [0, 0, m * g * l / 2]
    b = Bodies()
    ground = RBody!(b, 1, ones(3))

    p0 = normalize([1.0, 2, 3, 0])
    horizontal_beam = RBody!(b, 
        m,
        (m * l ^ 2) * [1 / 1200, 1 / 12, 1 / 12],
        [l / 2; 0.0; 0.0; p0],
    )

    s_rot = [0.0, 0, 0]
    rot_axis = [0.0, 0, 1]
    joints = [
        JointSimple(ground),  # Fix the ground
        JointPoint(ground, horizontal_beam, s_rot),
        JointPerpend1(ground, horizontal_beam, rot_axis, [1, 0, 0]),
        JointPerpend1(ground, horizontal_beam, rot_axis, [0, 1, 0]),
    ]

    # Add equivalent torque
    frc = ForceTorque(horizontal_beam, torque)

    sys = Mbs(b, joints, [frc], gv=[0, -g, 0])

    osys = OdeMbs(sys)
    
    # yp0 = sys(0.0, y0)
    y0 = initial_position(sys)
    dy = zeros(sys.ny)
    ode!(dy, y0, osys, 0.0)
    # acceleration and velocity should be zero
    @test norm(dy) < 1e-13
end