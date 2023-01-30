abstract type AbstractFrame end

mutable struct PointFrame <: AbstractFrame
    origin :: SVector{3}
end

mutable struct EulerParameterFrame <: AbstractFrame
    origin :: SVector{3}
    euler_parameters :: SVector{4}
end