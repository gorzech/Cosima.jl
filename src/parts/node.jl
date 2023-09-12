abstract type Node end

struct RBodyNode <: Node
    q0::SVector{7,Float64}
    h0::SVector{6,Float64}
    q::MVector{7,Float64}
    h::MVector{6,Float64}
    qi::UnitRange{Int}
    hi::UnitRange{Int}
end

nq(::RBodyNode) = 7
nh(::RBodyNode) = 6

rot(n::RBodyNode) = rot(n.q[4:7])

function point_global_to_local(n::RBodyNode, point)
    r = n.q[1:3]
    return global_to_local(n, point - r)
end

function global_to_local(n::RBodyNode, vector)
    A = rot(n)
    return Vector(A' * vector)
end