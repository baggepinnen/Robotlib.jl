using Robotlib
using Robotlib: Rt2T
function leds_to_tf(L::AbstractMatrix)
    size(L) == (3,3) || throw(ArgumentError("Expected a size 3×3 matrix"))
      m   = mean(L, dims=2)[:]
    # v1  = normalize(L[:,2]-L[:,1])
    # v2  = normalize(L[:,3]-L[:,1])
      v1  = normalize(L[:,1]-m)
      v2  = normalize(L[:,2]-m)
      v30 = normalize(L[:,2]-m)
      v3  = normalize(cross(v1, v2))
    v3'v30 < 0 && (v3 = -v3)
    v2 = normalize(cross(v3, v1))
    Rt2T([v1 v2 v3], m)
end

function leds_to_tf(L::AbstractVector{<:AbstractMatrix})
    length(L) == 3 || throw(ArgumentError("Expected three matrices"))
    N = size(L[1],1)
    T = Matrix{Float64}[]
    for i = 1:N
        D = [L[1][i,:] L[3][i,:] L[3][i,:]]
        all(isfinite, D) || continue
        push!(T, leds_to_tf(D))
    end
    T
end



"""
    T_RB_T, T_TF_TCP, normDist = calibAXYB(POSES, MEASURED, ortho)

This functions solves the problem AX = YB

- `POSES` is always the position of the tool frame from the robot FK
- `MEASURED` is the position of the object mounted on the tool frame, for
    instance if measuring with the Nikon it is the diod frame. If measuring
    with a laser scanner mounted on the TF, it is the position of the SCANNER
    in the measured object frame, not the object in the scanner frame!
- `ortho` is a flag indicating if found matrices are to be orthogonalized
"""
function calibAXYB(POSES, MEASURED, ortho = true)
    length(POSES) == length(MEASURED) ||
        throw(ArgumentError("POSES and MEASURED must have the same length"))
    T = promote_type(eltype(POSES[1]), eltype(MEASURED[1]))
    A = Matrix{T}[]
    B = Vector{T}[]
    I = I(3)
    Z = zeros(3)

    for i = 1:length(POSES)
        RA = POSES[i][1:3,1:3]
        BB = MEASURED[i][1:3,1:4]
        RB = BB[1:3,1:3]
        TB = BB[1:3,4]

        At = [RA Z  Z  Z    -I*RB[1,1] -I*RB[2,1] -I*RB[3,1]   Z
              Z  RA Z  Z    -I*RB[1,2] -I*RB[2,2] -I*RB[3,2]   Z
              Z  Z  RA Z    -I*RB[1,3] -I*RB[2,3] -I*RB[3,3]   Z
              Z  Z  Z  RA   -I*TB[1]   -I*TB[2]   -I*TB[3]    -I]

        b = [zeros(9,1); -POSES[i][1:3,4]]
        push!(A,At)
        push!(B,b)
    end
    A = reduce(vcat, A)
    B = reduce(vcat, B)

    w = A\B
    T_RB_T        = reshape(w[13:24],3,4)
    T_RB_T[4,:]   = [0 0 0 1]
    T_TF_TCP      = reshape(w[1:12],3,4)
    T_TF_TCP[4,:] = [0 0 0 1]
    if ortho
        T_RB_T   = Robotlib.orthogonalized(T_RB_T)
        T_TF_TCP = Robotlib.orthogonalized(T_TF_TCP)
        # Second optimization
        B2 = B-A[:,1:9]   * reshape(T_TF_TCP[1:3,1:3],9,1) - A[:,13:21] * reshape(T_RB_T[1:3,1:3]  ,9,1)
        A2 = A[:,[10:12 22:24]]
        w2 = A2\B2
        T_TF_TCP[1:3,4] = w2[1:3]
        T_RB_T[1:3,4]   = w2[4:6]
    end

    # Eval
    for i = 1:length(POSES)
        YN[i] = Robotlib.orthogonalized(T_RB_T*MEASURED[i])
        MX[i] = Robotlib.orthogonalized(POSES[i]*T_TF_TCP)
        dist  = (YN[i][1:3,4] - MX[i][1:3,4])'
        normDist[i] = norm(dist)
    end
    T_RB_T, T_TF_TCP, normDist
end