abstract type Joint end

struct JointPoint{T1<:Body,T2<:Body} <: Joint
    v_i::Vector{Float64}
    v_j::Vector{Float64}
    body_i::T1
    body_j::T2
end

nc(::JointPoint) = 3

function JointPoint(body_i, body_j, location)
    v_i = point_global_to_local(body_i, location)
    v_j = point_global_to_local(body_j, location)
    JointPoint(v_i, v_j, body_i, body_j)
end

function point_body_constr(b::RBody, v::Vector{Float64})
    A_i = rot(b.q[4:7])
    om_s = skew(b.h[4:6])
    s_i = v[1]
    c = b.q[1:3] + A_i * s_i
    c_q = [Matrix(I, 3, 3) -A_i * skew(s_i)]
    g = -A_i * om_s * om_s * s_i
    c_p = c_q * b.h
    return (c, c_q, c_p, g)
end

function constraints!(c, c_q, c_p, g, j::JointPoint)
    c_i, c_q_i, c_p_i, g_i = point_body_constr(j.body_i, j.v_i)
    c_j, c_q_j, c_p_j, g_j = point_body_constr(j.body_j, j.v_j)
    c .= c_i - c_j
    c_p .= c_p_i - c_p_j
    g .= g_i - g_j
    c_q[:, j.body_i.hi] = c_q_i
    c_q[:, j.body_j.hi] = -c_q_j
end

struct JointSimple{T<:Body} <: Joint
    body::T
    qi::Vector{Int}
    hi::Vector{Int}
    G0_2::Matrix{Float64}
end

nc(j::JointSimple) = length(j.qi)

function JointSimple(body, position_idx = [1, 2, 3], fix_rotation = true)
    @assert length(position_idx) <= 3
    @assert all(position_idx .>= 1)
    @assert all(position_idx .<= 3)
    if fix_rotation
        G0_2 = 0.5 * gep(body.q0[4:7])
        G0_2 = Matrix(G0_2[:, 2:end]') # To remove adjoint property
        qi = [position_idx; 5; 6; 7]
        hi = [position_idx; 4; 5; 6]
    else
        G0_2 = zeros(0, 0)
        qi = hi = position_idx
    end

    JointSimple(body, qi, hi, G0_2)
end

function constraints!(c, c_q, c_p, _, j::JointSimple)
    b = j.body
    c .= b.q[j.qi] .- b.q0[j.qi]
    # g is equal to zero
    if isempty(j.G0_2) # only tranlational part
        ncr = nc(j)
        c_q[:, b.hi[j.hi]] .= Matrix(I, ncr, ncr)
        c_p .= b.h[j.hi]
    else
        ncr = nc(j) - 3
        c_q[1:ncr, b.hi[j.hi[1:ncr]]] .= Matrix(I, ncr, ncr)
        c_q[ncr+1:end, b.hi[j.hi[ncr+1:end]]] .= j.G0_2
        c_p[1:ncr] .= b.h[j.hi[1:ncr]]
        c_p[ncr+1:end] .= j.G0_2 * b.h[j.hi[ncr+1:end]]
    end
end

struct JointPerpend1{T1<:Body,T2<:Body} <: Joint
    body_i::T1
    body_j::T2
    v_i::SVector{3,Float64}
    v_j::SVector{3,Float64}
end

nc(::JointPerpend1) = 1

function JointPerpend1(body_i, body_j, v_i, v_j, location = zeros(3))
    v_i = normalize(v_i)
    v_j = normalize(v_j)
    if v_i' * v_j > 1e-5
        throw(
            ArgumentError(
                "Vectors seems not to be perpendicular in JointPerpend1 constructor.",
            ),
        )
    end

    JointPerpend1(
        body_i,
        body_j,
        rot(body_i.q0[4:7])' * v_i, # global to local
        rot(body_j.q0[4:7])' * v_j,
    )
end

function body_vector_deriv(b::T, v_local) where {T<:RBody}
    om = b.h[4:6]
    om_s = skew(om)
    Ai = rot(b.q[4:7])
    vi = Ai * v_local
    cq_h = -Ai * skew(v_local)
    vi_p = cq_h * om
    vi_b_g = Ai * (om_s * (om_s * v_local))

    return (vi, cq_h, vi_p, vi_b_g)
end

function constraints!(c, c_q, c_p, g, j::JointPerpend1)
    v_i, cq_i, v_i_p, v_i_b_g = body_vector_deriv(j.body_i, j.v_i)
    v_j, cq_j, v_j_p, v_j_b_g = body_vector_deriv(j.body_j, j.v_j)
    c .= v_i' * v_j
    c_p .= v_i' * v_j_p + v_j' * v_i_p
    g .= -v_i' * v_j_b_g - v_j' * v_i_b_g - 2 * (v_i_p' * v_j_p)

    c_q[:, j.body_i.hi[4:end]] .= v_j' * cq_i
    c_q[:, j.body_j.hi[4:end]] .= v_i' * cq_j
    nothing
end