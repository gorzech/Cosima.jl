abstract type Joint end

struct JointPoint{N1<:Node,N2<:Node} <: Joint
    v_i::Vector{Float64}
    v_j::Vector{Float64}
    node_i::N1
    node_j::N2
end

nc(::JointPoint) = 3

JointPoint(body_i::Body, body_j::Body, location) = JointPoint(body_i.node, body_j.node, location)

function JointPoint(node_i::Node, node_j::Node, location)
    v_i = point_global_to_local(node_i, location)
    v_j = point_global_to_local(node_j, location)
    JointPoint(v_i, v_j, node_i, node_j)
end

function point_body_constr(n::RBodyNode, s_i::Vector{Float64})
    A_i = rot(n)
    om_s = skew(n.h[4:6])
    c = n.q[1:3] + A_i * s_i
    c_q = [Matrix(I, 3, 3) -A_i * skew(s_i)]
    g = -A_i * om_s * om_s * s_i
    c_p = c_q * n.h
    return (c, c_q, c_p, g)
end

function constraints!(c, c_q, c_p, g, j::JointPoint)
    c_i, c_q_i, c_p_i, g_i = point_body_constr(j.node_i, j.v_i)
    c_j, c_q_j, c_p_j, g_j = point_body_constr(j.node_j, j.v_j)
    c .= c_i - c_j
    c_p .= c_p_i - c_p_j
    g .= g_i - g_j
    c_q[:, j.node_i.hi] = c_q_i
    c_q[:, j.node_j.hi] = -c_q_j
end

struct JointSimple{N<:Node} <: Joint
    node::N
    qi::Vector{Int}
    hi::Vector{Int}
    G0_2::Matrix{Float64}
end

nc(j::JointSimple) = length(j.qi)

function JointSimple(body, position_idx=[1, 2, 3], fix_rotation=true)
    @assert length(position_idx) <= 3
    @assert all(position_idx .>= 1)
    @assert all(position_idx .<= 3)
    if fix_rotation
        G0_2 = 0.5 * gep(body.node.q0[4:7])
        G0_2 = Matrix(G0_2[:, 2:end]') # To remove adjoint property
        qi = [position_idx; 5; 6; 7]
        hi = [position_idx; 4; 5; 6]
    else
        G0_2 = zeros(0, 0)
        qi = hi = position_idx
    end

    JointSimple(body.node, qi, hi, G0_2)
end

function constraints!(c, c_q, c_p, _, j::JointSimple)
    n = j.node
    c .= n.q[j.qi] .- n.q0[j.qi]
    # g is equal to zero
    if isempty(j.G0_2) # only tranlational part
        ncr = nc(j)
        c_q[:, n.hi[j.hi]] .= Matrix(I, ncr, ncr)
        c_p .= n.h[j.hi]
    else
        ncr = nc(j) - 3
        c_q[1:ncr, n.hi[j.hi[1:ncr]]] .= Matrix(I, ncr, ncr)
        c_q[ncr+1:end, n.hi[j.hi[ncr+1:end]]] .= j.G0_2
        c_p[1:ncr] .= n.h[j.hi[1:ncr]]
        c_p[ncr+1:end] .= j.G0_2 * n.h[j.hi[ncr+1:end]]
    end
end

struct JointPerpend1{N1<:Node,N2<:Node} <: Joint
    node_i::N1
    node_j::N2
    v_i::SVector{3,Float64}
    v_j::SVector{3,Float64}

    function JointPerpend1(node_i::Node, node_j::Node, v_i, v_j)
        v_i = normalize(v_i)
        v_j = normalize(v_j)
        if v_i' * v_j > 1e-5
            throw(
                ArgumentError(
                    "Vectors seems not to be perpendicular in JointPerpend1 constructor.",
                ),
            )
        end

        new{typeof(node_i),typeof(node_j)}(
            node_i,
            node_j,
            global_to_local(node_i, v_i),
            global_to_local(node_j, v_j),
        )
    end
end

JointPerpend1(body_i::Body, body_j::Body, v_i, v_j) = JointPerpend1(body_i.node, body_j.node, v_i, v_j)

nc(::JointPerpend1) = 1

function body_vector_deriv(n::RBodyNode, v_local)
    om = n.h[4:6]
    om_s = skew(om)
    Ai = rot(n)
    vi = Ai * v_local
    cq_h = -Ai * skew(v_local)
    vi_p = cq_h * om
    vi_b_g = Ai * (om_s * (om_s * v_local))

    return (vi, cq_h, vi_p, vi_b_g)
end

function constraints!(c, c_q, c_p, g, j::JointPerpend1)
    v_i, cq_i, v_i_p, v_i_b_g = body_vector_deriv(j.node_i, j.v_i)
    v_j, cq_j, v_j_p, v_j_b_g = body_vector_deriv(j.node_j, j.v_j)
    c .= v_i' * v_j
    c_p .= v_i' * v_j_p + v_j' * v_i_p
    g .= -v_i' * v_j_b_g - v_j' * v_i_b_g - 2 * (v_i_p' * v_j_p)

    c_q[:, j.node_i.hi[4:end]] .= v_j' * cq_i
    c_q[:, j.node_j.hi[4:end]] .= v_i' * cq_j
    nothing
end