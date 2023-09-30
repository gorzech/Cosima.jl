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

function addbody!(sys::Mbs, body::RBody)
    append!(sys.bodies.r_bodies, (body,))
    sys.nq += nq(body)
    sys.nh += nh(body)
    sys.ny = sys.nq + sys.nh
end
