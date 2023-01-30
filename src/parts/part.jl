abstract type AbstractPart end

struct Rigid <: AbstractPart
    reference_frame :: AbstractFrame
    mass :: Float64
end

