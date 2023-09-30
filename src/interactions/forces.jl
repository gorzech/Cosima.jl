abstract type Force end

struct ForceTorque{T<:Body} <: Force
    body::T
    psi::Matrix{Float64}
    torque::Function
end

function ForceTorque(body::RBody, torque::Function)
    @assert size(torque(0.0)) == (3,)
    ForceTorque(body, zeros(0, 0), torque)
end

function force!(Q, t, f::ForceTorque{RBody})
    T = f.torque(t)
    pp_i = f.body.node.q[4:7]
    A_i = rot(pp_i)
    T_i = A_i' * T 

    Q[f.body.node.hi[4:end]] .+= T_i
end

function force(t, f::ForceTorque{<:Body})
    Q = empty(nh(f.body))
    force!(Q, t, f)
    return Q
end