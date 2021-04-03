using Robotlib
using Robotlib: Rt2T, T2R, T2t
function leds_to_tf(L::AbstractMatrix, usemean = false)
    size(L) == (3, 3) || throw(ArgumentError("Expected a size 3×3 matrix"))
    
    if usemean
        m = mean(L, dims = 2)[:] 
        v1 = normalize(L[:,1] - m)
        v2 = normalize(L[:,2] - m)
    else
        m = L[:,1]
        v1  = normalize(L[:,2]-m)
        v2  = normalize(L[:,3]-m)
    end
    v3 = normalize(cross(v1, v2))
    v2 = normalize(cross(v3, v1))
    R = [v1 v2 v3]
    det(R) < 0 && (R = [v1 v2 -v3])
    Rt2T(R, m)
end

function leds_to_tf(L::AbstractVector{<:AbstractMatrix})
    length(L) == 3 || throw(ArgumentError("Expected three matrices"))
    NL, N = size(L[1])
    NL > N && @warn "More LEDs than number of data points, are you sure you should not transpose the L matrices?"
    T = Matrix{Float64}[]
    for i = 1:N
        D = [L[1][:, i] L[2][:, i] L[3][:, i]]
        all(isfinite, D) || continue
        push!(T, leds_to_tf(D))
    end
    T
end



"""
    T_RB_T, T_TF_TCP, normDist = calibAXYB(POSES, MEASURED; ortho=true, estimator = \\, refine = true)

This functions solves the problem AX = YB

- `POSES` is always the position of the tool frame from the robot FK
- `MEASURED` is the position of the object mounted on the tool frame, for
    instance if measuring with the Nikon it is the diod frame. If measuring
    with a laser scanner mounted on the TF, it is the position of the SCANNER
    in the measured object frame, not the object in the scanner frame!
- `ortho` is a flag indicating if found matrices are to be orthogonalized
- `refine`: run a second, local optimization routine to improve the accuracy (recommended)
"""
function calibAXYB(POSES, MEASURED; ortho = true, estimator = \, refine = true)
    length(POSES) == length(MEASURED) ||
        throw(ArgumentError("POSES and MEASURED must have the same length"))
    T = promote_type(eltype(POSES[1]), eltype(MEASURED[1]))
    A = Matrix{T}[]
    B = Vector{T}[]
    I = Robotlib.I3
    Z = zeros(3, 3)

    for i = 1:length(POSES)
        RA = T2R(POSES[i])
        BB = MEASURED[i][1:3, 1:4]
        RB = T2R(BB)
        TB = T2t(BB)

        At = [
            RA Z Z Z -I*RB[1, 1] -I*RB[2, 1] -I*RB[3, 1] Z
            Z RA Z Z -I*RB[1, 2] -I*RB[2, 2] -I*RB[3, 2] Z
            Z Z RA Z -I*RB[1, 3] -I*RB[2, 3] -I*RB[3, 3] Z
            Z Z Z RA -I*TB[1] -I*TB[2] -I*TB[3] -I
        ]

        b = [zeros(9); -T2t(POSES[i])]
        push!(A, At)
        push!(B, b)
    end
    A = reduce(vcat, A)
    B = reduce(vcat, B)

    w = estimator(A, B)
    T_RB_T = reshape(w[13:24], 3, 4)
    T_RB_T = [T_RB_T; [0 0 0 1]]
    T_TF_TCP = reshape(w[1:12], 3, 4)
    T_TF_TCP = [T_TF_TCP; [0 0 0 1]]
    if ortho
        T_RB_T = Robotlib.orthonormal(T_RB_T)
        T_TF_TCP = Robotlib.orthonormal(T_TF_TCP)
        # Second optimization
        B2 = B - A[:, 1:9] * vec(T_TF_TCP[1:3, 1:3]) - A[:, 13:21] * vec(T_RB_T[1:3, 1:3])
        A2 = A[:, [10:12; 22:24]]
        w2 = estimator(A2, B2)
        T_TF_TCP[1:3, 4] = w2[1:3]
        T_RB_T[1:3, 4] = w2[4:6]
    end

    # Eval
    YN = similar(MEASURED)
    MX = similar(MEASURED)
    normDist = zeros(length(MEASURED))
    function cost(T_RB_T, T_TF_TCP)
        for i = 1:length(POSES)
            YN[i] = Robotlib.orthonormal(T_RB_T * MEASURED[i])
            MX[i] = Robotlib.orthonormal(POSES[i] * T_TF_TCP)
            dist = (YN[i][1:3, 4] - MX[i][1:3, 4])'
            normDist[i] = norm(dist)
        end
        sum(normDist)
    end
    cost(T_RB_T, T_TF_TCP)
    if refine
        p0 = zeros(12)
        adapter = function(p)
            T_RB_Ti = T_RB_T*Rt2T(rpy2R(p[4:6]), p[1:3])
            T_TF_TCPi = T_TF_TCP*Rt2T(rpy2R(p[10:12]), p[7:9])
            T_RB_Ti, T_TF_TCPi
        end
        cost(p) = cost(adapter(p)...)
        res = Optim.optimize(
            cost,
            p0,
            BFGS(),
            Optim.Options(
                store_trace       = true,
                show_trace        = true,
                show_every        = 2,
                iterations        = 100,
                allow_f_increases = false,
                time_limit        = 100,
                g_tol             = 1e-12,
            ),
        )
        T_RB_T, T_TF_TCP = adapter(res.minimizer)
        cost(T_RB_T, T_TF_TCP)
    end

    T_RB_T, T_TF_TCP, normDist
end

random_SE3(σt=1) = Rt2T(orthonormal(randn(3,3)), σt*randn(3))

"""
    T_RB_T, T_TF_S, POSES, MEASURED, LEDs = simulateCalibration_AXYB(N; σt=0.1, σT=1, σy=0)

σt denots the std of the translation of the sensor frame. σT is the std for the translation of the camera to the robot. σy is measurement noise
"""
function simulateCalibration_AXYB(N=100; σt=0.1, σT=1, σy=0)
    POSES = [random_SE3(0.5*σT) for _ = 1:N]
    T_RB_T = random_SE3(σT)
    T_TF_S = random_SE3(σt)
    MEASURED = [trinv(T_RB_T)*P*T_TF_S for P in POSES]
    LEDs = map(1:3) do l
        L = map(1:N) do i
            t = T2t(MEASURED[i])
            R = T2R(MEASURED[i])
            t .+ (l == 1 ? 0 : σt*R[:,l-1]) # origin + one of the basis vectors
        end
        reduce(hcat, L) .+ σy .* randn.()
    end
    T_RB_T, T_TF_S, POSES, MEASURED, LEDs
end
