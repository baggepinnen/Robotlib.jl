"""
    Rf, m = calib_force(POSES, F, m0::Real = 0.3; offset=true)
    Rf, m = calib_force(POSES, F, g::Vector; offset=true)

`POSES` ∈ ℜ(4,4,N) are the positions of the tool frame from the robot forward kinematics, 4x4 transformation matrices
`F` ∈ ℜ(N,3) vector of forces (accepts ℜ(Nx6) matrix with torques also)
# Usage
```
Rf*force[i,1:3] + forceoffs = POSES[1:3,1:3,i]'*[0, 0, mf*-9.82]
```
If `m0` is supplied, it's assumed that the gravity vector is [0,0,-g], or in words, that the gravity is acting along the negative z axis. A custom gravity vector `g` can also be supplied.

This function can also be used to calibrate an accelerometer.

[Bagge Carlson, F.](https://www.control.lth.se/staff/fredrik-bagge-carlson/), ["Machine Learning and System Identification for Estimation in Physical Systems"](https://lup.lub.lu.se/search/publication/ffb8dc85-ce12-4f75-8f2b-0881e492f6c0) (PhD Thesis 2018).

See also `calib_force_iterative` and `calib_force_eigen`.
"""
function calib_force(POSES,F,m0::Real=0.3; kwargs...)
    g0 = [0,0,-m0*9.82]
    calib_force(POSES,F,g0; kwargs...)
end
function calib_force(POSES,F,g::AbstractVector; offset=true, verbose=true)

    N = size(POSES,3)

    local forceoffs
    I = Robotlib.I3
    A  = Array{eltype(F)}(undef, 3N, offset ? 12 : 9)
    B  = Array{eltype(F)}(undef, 3N)
    A2 = Array{eltype(F)}(undef, 3N, offset ? 4 : 1)
    B2 = Array{eltype(F)}(undef, 3N)
    At = offset ? (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I   -I] : (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I]

    for i = 1:N
        RA               = POSES[1:3,1:3,i]'
        b                = RA*g
        A[3(i-1)+1:3i,:] = At(F,i)
        B[3(i-1)+1:3i]   = b
    end
    verbose && @info("Condition number of Gram matrix: ", cond(A'A))
    w  = tls(A,B)
    Rf = reshape(w[1:9],3,3)
    verbose && @info("Determinant of Rf before orthonormalization: ", det(Rf), opnorm(Rf))
    if det(Rf) < 0
        verbose && @warn("det(Rf) < 0, left handed coordinate system? Trying to solve the problem with reversed gravity (-g)")
        return calib_force(POSES,F,-g; offset=offset, verbose=verbose)
    end
    orthonormal!(Rf)
    At = offset ? RA -> [RA*g/norm(g) -I] : RA -> RA*g/norm(g)

    for i = 1:N
        RA                = POSES[1:3,1:3,i]'
        b                 = Rf*F[i,1:3]
        A2[3(i-1)+1:3i,:] = At(RA)
        B2[3(i-1)+1:3i]   = b
    end
    w         = A2\B2
    m         = w[1]/9.82
    offset && (forceoffs = w[2:4])

    m < 0 && @error("Estimated mass is negative, left handed coordinate system for sensor?")
    if offset
        return Rf,m,forceoffs
    else
        return Rf, m
    end
end





"""
    Rf, g, m = calib_force_eigen(POSES, F, g)

`POSES` ∈ ℜ(4,4,N) or ℜ(3,3,N) are the orientations of the tool frame from the robot forward kinematics, 4×4 transformation matrices or 3×3 rotation matrices.
`F` ∈ ℜ(N,3) vector of forces (accepts ℜ(N×6) matrix with torques also)
´g´ is an initial guess for the gravity vector.
# Usage
```
Rf*force[i,1:3] = POSES[1:3,1:3,i]'*g
```
See also `calib_force_eigen`.
"""
function calib_force_iterative(POSES, F, g; trace = false)
    N = size(POSES, 3)
    I = Robotlib.I3
    A = Array{eltype(F)}(undef, 3N, 3)
    B = F'[:]
    local m, Rf
    trace && (Rg = [])
    for iter = 1:6
        Rf, m = calib_force(
            POSES,
            F,
            g;
            offset = false,
            verbose = false,
        )
        for i = 1:N
            A[3(i-1)+1:3i, :] = Rf'POSES[:, :, i]'
        end
        g = A \ B
        trace && push!(Rg, (Rf, g))
    end
    trace && (return Rf, g, m, Rg)
    Rf, g, m
end

function eigenR(K)
    v = real(eigen(K).vectors[:, end])
    R = reshape(v, 3, 3)
    det(R) < 0 && (R .*= -1)
    orthonormal(R)
end

"""
    calib_force_eigen(POSES, F)

`POSES` ∈ ℜ(4,4,N) or ℜ(3,3,N) are the orientations of the tool frame from the robot forward kinematics, 4×4 transformation matrices or 3×3 rotation matrices.
`F` ∈ ℜ(N,3) vector of forces (accepts ℜ(N×6) matrix with torques also)
# Usage
```
Rf*force[i,1:3] + forceoffs = POSES[1:3,1:3,i]'*[0, 0, mf*-9.82]
```
"""
function calib_force_eigen(POSES, forces)
    N = size(POSES, 3)
    D = reshape(POSES[1:3,1:3,:], 3, 3N)'
    F = kron(forces, I(3))
    K = F \ D * (D \ F)
    R̂ = eigenR(K)
    ĝ = D \ F * vec(R̂)
    R̂, ĝ
end

import Robotlib: skew

"""
    This function uses Cayley transform to solve for R. It requires a known mass so it is recommended to use `calib_force` instead.
"""
function calib_force2(POSES,F,m,g  = -9.82)
    N  = size(POSES,3)
    
    mg = m*g

    I  = Robotlib.I3
    A  = Array{eltype(F)}(undef, 3N, 3)
    B  = Array{eltype(F)}(undef, 3N)

    for i = 1:N
        RA               = POSES[3,1:3,i]
        A[3(i-1)+1:3i,:] = skew(F[i,1:3][:] + RA*mg)
        B[3(i-1)+1:3i]   = RA*mg - F[i,1:3]
    end
    @info("Condition number of Gram matrix: ", cond(A'A))
    w  = tls(A,B)
    S  = skew(w)
    Rf = (I+S)\(I-S)
end
