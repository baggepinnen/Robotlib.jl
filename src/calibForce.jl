"""
`calibForce(POSES, F, m0=0.3; offset=true)`
If result is bad, check if you send data in correct form;
`POSES` ∈ ℜ(4,4,N) is always the position of the tool frame from the robot FK,
4x4 transformation matrices
`F` ∈ ℜ(N,3) vector of forces (accepts ℜ(Nx6) matrix with torques also)
usage `Rf*force[i,1:3] + forceoffs = POSES[1:3,1:3,i]'*[0, 0, mf*-9.82]`
Bagge
"""
function calibForce(POSES,F,m0=0.3; offset=true)

    N = size(POSES,3)

    local forceoffs
    g = 9.82
    mg = m0*g # Initial guess is m0, the method accepts 5 orders of magnitude error at least

    I = eye(3)
    A  = Array{eltype(F)}(undef, 3N, offset ? 12 : 9)
    B  = Array{eltype(F)}(undef, 3N)
    A2 = Array{eltype(F)}(undef, 3N, offset ? 4 : 1)
    B2 = Array{eltype(F)}(undef, 3N)
    At = offset ? (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I   I] : (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I]

    for i = 1:N
        RA               = POSES[1:3,1:3,i]'
        b                = -RA[1:3,3]
        A[3(i-1)+1:3i,:] = At(F,i)
        B[3(i-1)+1:3i]   = b
    end
    @info("Condition number of Gram matrix: ", cond(A'A))
    w  = A\B.*mg
    Rf = reshape(w[1:9],3,3)
    @info("Determinant of Rf before orthonormalization: ", det(Rf)," (Should be close to 1, but don't be afraid if it's not, should however be positive!)")
    det(Rf) < 0 && @error("det(Rf) < 0, left handed coordinate system?")
    toOrthoNormal!(Rf)
    At = offset ? RA -> [RA[:,3] I] : RA -> RA[:,3]
    for ii = 1:2
        for i = 1:N
            RA                = POSES[1:3,1:3,i]'
            b                 = -Rf*F[i,1:3]
            A2[3(i-1)+1:3i,:] = At(RA)
            B2[3(i-1)+1:3i]   = b
        end
        w         = A2\B2
        mg        = w[1]
        offset && (forceoffs = w[2:4])
    end
    m = mg/g
    m < 0 && @error("Estimated mass is negative, left handed coordinate system for sensor?")
    if offset
        return Rf,m,forceoffs
    else
        return Rf,m
    end
end
