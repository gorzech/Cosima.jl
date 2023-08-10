skew(x) = SA[0.0 -x[3] x[2]; x[3] 0.0 -x[1]; -x[2] x[1] 0.0]

gep(p) = SA[
    -p[2] p[1] p[4] -p[3]
    -p[3] -p[4] p[1] p[2]
    -p[4] p[3] -p[2] p[1]
]

eep(p) = SA[
    -p[2] p[1] -p[4] p[3]
    -p[3] p[4] p[1] -p[2]
    -p[4] -p[3] p[2] p[1]
]

# ROT returns 3x3 rotation matrix for a vector of Euler parameters"""
rot(p) = eep(p) * gep(p)'