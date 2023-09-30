abstract type Body end

struct RBody <: Body
    node::RBodyNode
    mass::Float64
    Ic::SMatrix{3,3,Float64}
end

nq(b::RBody) = nq(b.node)
nh(b::RBody) = nh(b.node)

struct Bodies
    r_bodies::Vector{RBody}
end

function Bodies()
    Bodies(RBody[])
end

function update_body_coordinates!(bodies::Bodies, y, nq)
    function update_coords!(n)
        p = y[n.qi[4:7]]
        normalize!(p)
        n.q[1:3] .= y[n.qi[1:3]]
        n.q[4:7] .= p
        n.q[8:end] .= y[n.qi[8:end]]
        n.h .= y[nq.+n.hi]
    end
    for b in bodies.r_bodies
        update_coords!(b.node)
    end
end

function last_body_idx(bodies::Bodies)
    qi, hi = 0, 0
    if !isempty(bodies.r_bodies)
        qi = bodies.r_bodies[end].node.qi[end]
        hi = bodies.r_bodies[end].node.hi[end]
    end
    return (qi, hi)
end

function set_body_coordinates!(b::Body, q)
    b.node.q .= q
end

function set_body_coordinates!(b::Body, q, h)
    b.node.q .= q
    b.node.h .= h
end

function set_body_coordinates!(b::Bodies, q)
    for bi in b.r_bodies
        bi.node.q .= q[bi.node.qi]
    end
end

function set_body_coordinates!(b::Bodies, q, h)
    for bi in b.r_bodies
        bi.node.q .= q[bi.node.qi]
        bi.node.h .= h[bi.node.hi]
    end
end

function mass_upper!(M, b::RBody)
    for i = 1:3
        M[i, i] = b.mass
    end
    M[4:6, 4:6] .= b.Ic
    nothing
end

function mass(b::Body)
    M = zeros(nh(b), nh(b))
    mass_upper!(M, b)
    return Symmetric(M, :U)
end

function grav(b::RBody, gv)
    return [gv .* b.mass; 0; 0; 0]
end

function quadratic_inertia(b::RBody)
    return [0.0; 0; 0; cross(b.node.h[4:end], b.Ic * b.node.h[4:end])]
end
