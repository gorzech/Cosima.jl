mutable struct EulerParameterFrame
    origin :: SVector{3, Float64}
    euler_parameters :: SVector{4, Float64}
end