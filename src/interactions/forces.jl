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

function ForceTorque!(sys::Mbs, body::RBody, torque::Function)
    ft = ForceTorque(body, zeros(0, 0), torque)
    append!(sys.forces, (ft,))
    return ft
end

function force!(Q, t, f::ForceTorque{RBody})
    T = f.torque(t)
    pp_i = f.body.q[4:7]
    A_i = rot(pp_i)
    T_i = A_i' * T 

    Q[f.body.hi[4:end]] .+= T_i
end

function force(t, f::ForceTorque{<:Body})
    Q = empty(f.body.nh)
    force!(Q, t, f)
    return Q
end