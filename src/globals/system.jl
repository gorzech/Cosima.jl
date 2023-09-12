struct Mbs
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

function Mbs(
    bodies::Bodies,
    joints::Vector{<:Joint},
    forces::Vector{<:Force};
    gv = SA[0.0, 0.0, 0.0],
    use_quadratic = true,
    baumg_params = (20.0, 20.0),
)
    nq, nh = last_body_idx(bodies)
    if norm(gv) > 0.0 || use_quadratic
        forces = convert(Vector{Force}, forces)
        if norm(gv) > 0.0
            append!(forces, (MbsGravForce(gv, bodies),))
        end
        if use_quadratic
            append!(forces, (MbsQuadraticForce(bodies),))
        end
    end
    nconstr = 0
    for j in joints
        nconstr += nc(j)
    end

    Mbs(
        bodies,
        joints,
        forces,
        (2 * baumg_params[1], baumg_params[2]^2),
        nq,
        nh,
        nq + nh,
        nconstr,
    )
end
