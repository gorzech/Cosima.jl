@enum axis x = 1 y = 2 z = 3

skew(x) = SA[0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]

gep(p) = SA[
    -p[2] p[1] p[4] -p[3]
    -p[3] -p[4] p[1] p[2]
    -p[4] p[3] -p[2] p[1]
]

eep(p) = SA[
    -p[2] p[1] -p[4] p[3]
    -p[3] p[4] p[1] -p[2]
    -p[4] -p[3] p[2] p[1]
]

# ROT returns 3x3 rotation matrix for a vector of Euler parameters"""
rot(p) = eep(p) * gep(p)'

function rot_axis(angle, around::axis = x)
    """ROTAXIS rotate by an angle (in RAD) around axis ('x', 'y', or 'z')"""
    ca = cos(angle)
    sa = sin(angle)
    res = zeros(3, 3)

    if around === x
        res[1, 1] = 1
        res[2, 2] = ca
        res[3, 3] = ca
        res[2, 3] = -sa
        res[3, 2] = sa
    elseif around === y
        res[2, 2] = 1
        res[1, 1] = ca
        res[3, 3] = ca
        res[1, 3] = sa
        res[3, 1] = -sa
    elseif around === z
        res[3, 3] = 1
        res[1, 1] = ca
        res[1, 2] = -sa
        res[2, 1] = sa
        res[2, 2] = ca
    else
        throw(ArgumentError("Wrong axis specification in call to rot_axis"))
    end
    return res
end

function q_axis(angle::Float64, around = axis::x)
    """Q_AXIS Quaternions that rotate by alpha (in RAD) around axis ('x', 'y', or 'z')"""
    p = [cos(angle / 2), 0, 0, 0]
    sa2 = sin(angle / 2)
    if around === x
        p[2] = sa2
    elseif around === y
        p[3] = sa2
    elseif around === z
        p[4] = sa2
    else
        throw(ArgumentError("Wrong axis specification in call to q_axis"))
    end
    return p
end

function q_mul(p1, p2)
    """Multiplication of quaternions. For rotational quat. rotation p2
    followed by p1: A(p1) * A(p2)"""
    [
        p1[1] * p2[1] - p1[2] * p2[2] - p1[3] * p2[3] - p1[4] * p2[4],
        p1[1] * p2[2] + p1[2] * p2[1] + p1[3] * p2[4] - p1[4] * p2[3],
        p1[1] * p2[3] - p1[2] * p2[4] + p1[3] * p2[1] + p1[4] * p2[2],
        p1[1] * p2[4] + p1[2] * p2[3] - p1[3] * p2[2] + p1[4] * p2[1],
    ]
end

function approximate_jacobian(vector_fun, x0, ε = 1e-8)
    # Approximate Jacobian matrix based on the finite differences
    # eps and delta influenced by Numerical Recipes by Press et. al
    # ε = 1e-8 # Approximately sqrt of machine precision
    x_Δ = copy(x0)
    f0 = vector_fun(x_Δ)
    res = zeros(length(f0), length(x0))
    for i = 1:length(x0)
        temp = x_Δ[i]
        Δ = ε * abs(temp)
        if Δ ≤ 0.0
            Δ = ε
        end
        x_Δ[i] = temp + Δ
        f_Δ = vector_fun(x_Δ)
        res[:, i] = (f_Δ - f0) / Δ # Forward difference formula
        x_Δ[i] = temp
    end
    return res
end