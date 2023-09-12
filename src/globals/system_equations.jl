function mass_upper!(M, mbs::Mbs)
    for b in mbs.bodies.r_bodies
        @views mass_upper!(M[b.node.hi, b.node.hi], b)
    end
end

function mass(mbs::Mbs)
    M = zeros(mbs.nh, mbs.nh)
    mass_upper!(M, mbs)
    return Symmetric(M, :U)
end

function force!(Q, _, fg::MbsGravForce)
    for b in fg.bodies.r_bodies
        Q[b.node.hi] .+= grav(b, fg.gv)
    end
end

function force!(Q, _, fqi::MbsQuadraticForce)
    for b in fqi.bodies.r_bodies
        Q[b.node.hi] .-= quadratic_inertia(b)
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
        y0[b.node.qi] .= b.node.q0
        y0[mbs.nq.+b.node.hi] .= b.node.h0
    end
    return y0
end

function compute_q_dot(qp, y, sys::Mbs)
    nq = sys.nq
    function body_q_dot(b)
        qp[b.node.qi[1:3]] .= y[nq.+b.node.hi[1:3]]
        om = y[nq.+b.node.hi[4:6]]
        p = b.node.q[4:7]
        L = gep(p)
        qp[b.node.qi[4:7]] .= 0.5 .* (L' * om) # p_dot
        qp[b.node.qi[8:end]] .= y[nq.+b.node.hi[7:end]]
    end
    for b in sys.bodies.r_bodies
        body_q_dot(b)
    end
end