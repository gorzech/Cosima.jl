mutable struct Mbs
    bodies::Bodies
    joints::Vector{<:Joint}
    forces::Vector{<:Force}
    baumg_params::Tuple{Float64,Float64}
    nq::Int
    nh::Int
    ny::Int
    nconstr::Int
end

struct MbsGravForce <: Force
    gv::SVector{3,Float64}
    bodies::Bodies
end

struct MbsQuadraticForce <: Force
    bodies::Bodies
end

function Mbs(;
        gv = SA[0.0, 0.0, 0.0],
        use_quadratic = true,
        baumg_params = (20.0, 20.0),
    )
    forces = Force[]
    if norm(gv) > 0.0
        append!(forces, (MbsGravForce(gv, bodies),))
    end
    if use_quadratic
        append!(forces, (MbsQuadraticForce(bodies),))
    end
    Mbs(
        Bodies(),
        Joint[],
        forces,
        (2 * baumg_params[1], baumg_params[2]^2),
        0,
        0,
        0,
        0,
    )
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
    node = RBodyNode(q0, h0, copy(q0), copy(h0), qi+1:qi+7, hi+1:hi+6)
    b = RBody(node, mass, Ic)
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

function addbody!(sys::Mbs, body::RBody)
    append!(sys.bodies.r_bodies, (body,))
    sys.nq += nq(body)
    sys.nh += nh(body)
    sys.ny = sys.nq + sys.nh
end

function ForceTorque!(sys::Mbs, body::RBody, torque::Function)
    ft = ForceTorque(body, zeros(0, 0), torque)
    append!(sys.forces, (ft,))
    return ft
end
