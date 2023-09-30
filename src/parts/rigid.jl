abstract type Body end

struct RBody <: Body
    nq::Int
    nh::Int
    q0::SVector{7,Float64}
    h0::SVector{6,Float64}
    q::MVector{7,Float64}
    h::MVector{6,Float64}
    qi::UnitRange{Int}
    hi::UnitRange{Int}
    mass::Float64
    Ic::SMatrix{3,3,Float64}
end

struct Bodies
    r_bodies::Vector{RBody}
end

function Bodies()
    Bodies(RBody[])
end

function last_body_idx(bodies::Bodies)
    qi, hi = 0, 0
    if !isempty(bodies.r_bodies)
        qi = bodies.r_bodies[end].qi[end]
        hi = bodies.r_bodies[end].hi[end]
    end
    return (qi, hi)
end

function set_body_coordinates!(b::Body, q)
    b.q .= q
end

function set_body_coordinates!(b::Body, q, h)
    b.q .= q
    b.h .= h
end

function set_body_coordinates!(b::Bodies, q)
    for bi in b.r_bodies
        bi.q .= q[bi.qi]
    end
end

function set_body_coordinates!(b::Bodies, q, h)
    for bi in b.r_bodies
        bi.q .= q[bi.qi]
        bi.h .= h[bi.hi]
    end
end


function RBody!(
    sys::Mbs,
    mass,
    Ic::Union{Matrix{Float64},SMatrix{3,3,Float64}},
    q0,
    h0,
)
    @assert mass > 0 "Mass must be positive"
    @assert size(Ic) == (3, 3) "Wrong size of an inertia matrix"
    @assert all(diag(Ic) .> 0) "Inertia diagonal terms must be positive"
    qi, hi = last_body_idx(sys.bodies)
    b = RBody(7, 6, q0, h0, copy(q0), copy(h0), qi+1:qi+7, hi+1:hi+6, mass, Ic)
    addbody!(sys, b)
    return b
end

function RBody!(sys::Mbs, mass, Ic::Union{Vector{Float64},SVector{3,Float64}}, q0, h0)
    RBody!(sys, mass, SMatrix{3,3}(diagm(Ic)), q0, h0)
end

function RBody!(sys::Mbs, mass, Ic, q0)
    RBody!(sys, mass, Ic, SVector{7}(q0), SVector{6}(zeros(6)))
end

function RBody!(sys::Mbs, mass, Ic)
    RBody!(sys, mass, Ic, SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], SVector{6}(zeros(6)))
end

function mass_upper!(M, b::RBody)
    for i = 1:3
        M[i, i] = b.mass
    end
    M[4:6, 4:6] .= b.Ic
    nothing
end

function mass(b::Body)
    M = zeros(b.nh, b.nh)
    mass_upper!(M, b)
    return Symmetric(M, :U)
end

function grav(b::RBody, gv)
    return [gv .* b.mass; 0; 0; 0]
end

function quadratic_inertia(b::RBody)
    return [0.0; 0; 0; cross(b.h[4:end], b.Ic * b.h[4:end])]
end

function point_global_to_local(b::RBody, point)
    r = b.q[1:3]
    A = rot(b.q[4:end])
    return Vector(A' * (point - r))
end

