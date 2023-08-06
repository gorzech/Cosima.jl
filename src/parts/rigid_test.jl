# Test version that follows the pattern where we have
# nodes - for coordinates delivery
# parts - for equation delivery
# frames - for interface delivery
# 
# Some rules
# each part has at least one node 
# frames can be attached to parts
# do we have always a CM frame? 
# how to apply gravity to a body - as we need an inertia info.

abstract type Node end
abstract type NodeFrame <: Node end

struct NodeQuaternions1 <: NodeFrame
    position_idx :: UnitRange{Int}
    velocity_idx :: UnitRange{Int}
end

NodeQuaternions = NodeQuaternions1

struct Rigid1
    nodecm :: NodeQuaternions
    mass :: Float64
    inertiacm :: SVector{6, Float64}
end

Rigid = Rigid1

mutable struct System1
    position :: Vector{float}
    velocity :: Vector{float}
    # nodes :: Vector{<:Node}
end

System = System1