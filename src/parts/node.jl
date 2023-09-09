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