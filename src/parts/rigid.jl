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


function RBody!(
    bodies::Bodies,
    mass,
    Ic::Union{Matrix{Float64},SMatrix{3,3,Float64}},
    q0,
    h0,
)
    @assert mass > 0 "Mass must be positive"
    @assert size(Ic) == (3, 3) "Wrong size of an inertia matrix"
    @assert all(diag(Ic) .> 0) "Inertia diagonal terms must be positive"
    qi, hi = last_body_idx(bodies)
    node = RBodyNode(q0, h0, copy(q0), copy(h0), qi+1:qi+7, hi+1:hi+6)
    b = RBody(node, mass, Ic)
    append!(bodies.r_bodies, (b,))
    return b
end

function RBody!(bodies::Bodies, mass, Ic::Union{Vector{Float64},SVector{3,Float64}}, q0, h0)
    RBody!(bodies, mass, SMatrix{3,3}(diagm(Ic)), q0, h0)
end

function RBody!(bodies::Bodies, mass, Ic, q0)
    RBody!(bodies, mass, Ic, SVector{7}(q0), SVector{6}(zeros(6)))
end

function RBody!(bodies::Bodies, mass, Ic)
    RBody!(bodies, mass, Ic, SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], SVector{6}(zeros(6)))
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
