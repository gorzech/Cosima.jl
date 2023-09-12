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

function acc_lambda_dae_index1!(p::OdeMbs, t)
    p.M .= 0.0 # clear the remaining due to factorization
    mass_upper!(p.M, p.mbs)
    force!(p.Rhs, t, p.mbs)
    if p.mbs.nconstr > 0
        constraints!(p.C, p.Cq, p.C′, p.G, p.mbs)
        p.M[1:p.mbs.nh, p.mbs.nh+1:p.mbs.nh+p.mbs.nconstr] .= p.Cq'
        p.Rhs[p.mbs.nh+1:end] .= p.G .- p.mbs.baumg_params[1] .* p.C′ .- p.mbs.baumg_params[2] .* p.C
    end
    bk = bunchkaufman!(Symmetric(p.M, :U))
    ldiv!(p.AccLambda, bk, p.Rhs)
end
