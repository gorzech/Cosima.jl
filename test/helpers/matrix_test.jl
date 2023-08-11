@testset "approximate_jacobian" begin
    """Basic test to evaluate approxJacobian"""

    function vector_test_function(x)
        x1, x2, x3 = x

        return [
            0.0
            146.0
            x1 + 4.50 * x2 * x2 + sqrt(x3) / 12.0
            13.0 * sin(4.0 * x1) + cos(3.0 * x2)
            1.0 / (x1 + x2 + x3)
            exp(2.0 * x1) + log(x2)
        ]
    end

    x0 = [1.23, 6.15, 9.0]
    jac = approximate_jacobian(vector_test_function, x0)

    row_4 = -1.0 / (x0[1] + x0[2] + x0[3])^2.0
    jac_exact = [
        0.0 0.0 0.0
        0.0 0.0 0.0
        1.0 9.0*x0[2] 1.0/24.0/sqrt(x0[3])
        52.0*cos(4.0 * x0[1]) -3.0*sin(3.0 * x0[2]) 0.0
        row_4 row_4 row_4
        2.0*exp(2.0 * x0[1]) 1.0/x0[2] 0.0
    ]

    @test jac ≈ jac_exact rtol = 1e-7 atol = 1e-5
end

@testset "q_mul_basic" begin
    p1 = normalize(rand(4))
    p2 = normalize(rand(4))

    pn = q_mul(p1, p2)

    # Now the same with rotation matrices
    A1 = rot(p1)
    A2 = rot(p2)
    An = rot(pn)

    A_expected = A1 * A2

    @test A_expected ≈ An rtol = 1e-14 atol = 1e-15
end
