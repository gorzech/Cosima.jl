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

function mass_upper!(M, mbs::Mbs)
    for b in mbs.bodies.r_bodies
        @views mass_upper!(M[b.hi, b.hi], b)
    end
end

function mass(mbs::Mbs)
    M = zeros(mbs.nh, mbs.nh)
    mass_upper!(M, mbs)
    return Symmetric(M, :U)
end

function force!(Q, _, fg::MbsGravForce)
    for b in fg.bodies.r_bodies
        Q[b.hi] .+= grav(b, fg.gv)
    end
end

function force!(Q, _, fqi::MbsQuadraticForce)
    for b in fqi.bodies.r_bodies
        Q[b.hi] .-= quadratic_inertia(b)
    end
end

function force!(Q, t, mbs::Mbs)
    Q .= 0
    for f in mbs.forces
        force!(Q, t, f)
    end
end

function constraints!(c, c_q, c_p, g, mbs::Mbs)
    c_idx = 0
    for j in mbs.joints
        jnc = nc(j)
        @views constraints!(
            c[c_idx+1:c_idx+jnc],
            c_q[c_idx+1:c_idx+jnc, :],
            c_p[c_idx+1:c_idx+jnc],
            g[c_idx+1:c_idx+jnc],
            j,
        )
        c_idx += jnc
    end
    nothing
end

function constraints(mbs::Mbs)
    c, cq, c′, g = zeros(mbs.nconstr),
    zeros(mbs.nconstr, mbs.nh),
    zeros(mbs.nconstr),
    zeros(mbs.nconstr)
    constraints!(c, cq, c′, g, mbs)
    return (c, cq, c′, g)
end

function force(t, mbs::Mbs)
    Q = zeros(mbs.nh)
    force!(Q, t, mbs)
    return Q
end

function initial_position(mbs::Mbs)
    y0 = zeros(mbs.ny)
    for b in mbs.bodies.r_bodies
        y0[b.qi] .= b.q0
        y0[mbs.nq.+b.hi] .= b.h0
    end
    return y0
end

function update_body_coordinates!(bodies::Bodies, y, nq)
    function update_coords!(b)
        p = y[b.qi[4:7]]
        normalize!(p)
        b.q[1:3] .= y[b.qi[1:3]]
        b.q[4:7] .= p
        b.q[8:end] .= y[b.qi[8:end]]
        b.h .= y[nq.+b.hi]
    end
    for b in bodies.r_bodies
        update_coords!(b)
    end
end

function compute_q_dot(qp, y, sys::Mbs)
    nq = sys.nq
    function body_q_dot(b)
        qp[b.qi[1:3]] .= y[nq.+b.hi[1:3]]
        om = y[nq.+b.hi[4:6]]
        p = b.q[4:7]
        L = gep(p)
        qp[b.qi[4:7]] .= 0.5 .* (L' * om) # p_dot
        qp[b.qi[8:end]] .= y[nq.+b.hi[7:end]]
    end
    for b in sys.bodies.r_bodies
        body_q_dot(b)
    end
end

struct OdeMbs
    mbs::Mbs
    M::Matrix{Float64}
    Rhs::Vector{Float64}
    C::Vector{Float64}
    Cq::Matrix{Float64}
    C′::Vector{Float64}
    G::Vector{Float64}
    AccLambda::Vector{Float64}
end

function OdeMbs(mbs::Mbs)
    OdeMbs(
        mbs,
        zeros(mbs.nh + mbs.nconstr, mbs.nh + mbs.nconstr),
        zeros(mbs.nh + mbs.nconstr),
        zeros(mbs.nconstr),
        zeros(mbs.nconstr, mbs.nh),
        zeros(mbs.nconstr),
        zeros(mbs.nconstr),
        zeros(mbs.nh + mbs.nconstr),
    )
end

function ode!(du, u, p::OdeMbs, t)
    update_body_coordinates!(p.mbs.bodies, u, p.mbs.nq)
    compute_q_dot(du, u, p.mbs)
    acc_lambda_dae_index1!(p, t)
    du[p.mbs.nq+1:end] .= p.AccLambda[1:p.mbs.nh]
    nothing
end
