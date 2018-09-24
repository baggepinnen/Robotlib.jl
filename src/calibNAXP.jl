import Robotlib: T2R, Rt2T
"""
    T_TF_S, meanRMS,RMS, norms = calibNAXP(points_S, lines_S, POSES, T_TF_S, planes::AbstractVector{Int},  iters = 50; doplot=false)

This functions implements the algorithm from the paper
"Six DOF eye-to-hand calibration from 2D measurements using planar constraints"
which solves the problem `Prb = n'A(X*Ps)`
If result is bad, check if you send data in correct form\n
`POSES ∈ R(4,4,N)` is always the position of the tool frame from the robot FK\n
`points_S ∈ R(3,N)` are measurements from a line laser scanner. Given
in the sensor frame. This script assumes that the laser plane is the
sensor XY plane with y-axis pointing outwards from the sensor.\n
`lines_S ∈ R(3,N)` are vectors corresponding to the direction of the measured laser line
`planes ∈ Z(N)` is a vector of indices corresponding to which plane a mesurement
comes from. It must be same length as `points_S` and enumerate the planes starting at 1. Example (N=15, number of planes = 3): [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]\n
`iters` determines the number of iterations\n
@inproceedings{carlson2015six,
  title={Six DOF eye-to-hand calibration from 2D measurements using planar constraints},
  author={Carlson, Fredrik Bagge and Johansson, Rolf and Robertsson, Anders},
  booktitle={Intelligent Robots and Systems (IROS), 2015 IEEE/RSJ International Conference on},
  pages={3628--3632},
  year={2015},
  organization={IEEE}
}
"""
function calibNAXP(points_S, lines_S, POSES, T_TF_S, planes::AbstractVector{Int},  iters::Integer = 50; doplot=false, trueT_TF_S=nothing)
    N_poses   = size(POSES,3)
    N_planes  = maximum(planes)
    RMScalibs = zeros(iters)
    ALLcalibs = zeros(iters, N_planes)
    norms     = trueT_TF_S == nothing ? nothing : zeros(iters)
    for c = 1:iters
        # Convert points to RB cordinate system
        points, lines = transform_data(POSES, points_S, lines_S, T_TF_S)
        normals, N_RB = findplanes(planes, points)

        # Estimation
        A,y = get_matrices(N_RB, POSES, points_S, lines_S)
        w   = c <= 5 ? A\y : tls(A,y)
        # w   = tls(A,y)
        R,t = w2Rt(w)

         # Re-estimation of translation
        w₂     = R[1:3,1:2][:]
        y₂     = y - A[:,1:6]*w₂
        t₂     = c <= 5 ? A[:,7:9]\y₂ : tls(A[:,7:9],y₂)
        T_TF_S = Rt2T(R,t₂)

        # Calculate iteration errors
        RMSi = map(1:N_planes) do j
            ind     = planes .== j
            ind     = findall(ind)
            sqrt(pointDiff(T_TF_S,POSES[:,:,ind],points_S[1:3,ind])[1])
        end
        RMScalibs[c]   = mean(RMSi)
        ALLcalibs[c,:] = RMSi
        trueT_TF_S == nothing || (norms[c] = norm(T_TF_S-trueT_TF_S))
        #  any(abs(RMSi) > 1e-1) && warn("Points does not seem to lie on a plane")
    end
    if doplot
        plotPlanes(normals)
        plotLines(points,lines)
    end
    T_TF_S, RMScalibs,ALLcalibs, norms
end

function w2Rt(w)
    H   = reshape(w,3,3)
    Rx  = H[1:3,1]
    Ry  = H[1:3,2]
    Rz  = cross(Rx,Ry) # These lines must be changed if laser plane is changed
    Rnn = [Rx Ry Rz]
    # @show Rnn
    R   = toOrthoNormal(Rnn)
    t   = H[1:3,3]
    R,t
end

function get_matrices(N_RB, POSES, points_S, lines_S)
    N_poses = size(POSES,3)
    A = zeros(2*N_poses,9)
    y = zeros(2*N_poses)
    for i = 1:N_poses
        Ni     = N_RB[:,i]
        Ra     = POSES[1:3,1:3,i]
        Ta     = POSES[1:3,4,i]
        Pi     = points_S[:,i]
        y[i]   = norm(Ni)^2-Ni'Ta
        A[i,:] = (reshape(repeat(Ni'Ra,3,1)',9,1).*[reshape(repeat(Pi[1:2],1,3)',6,1);1;1;1])' #TODO: rewrite
        # Sample one additional point along line
        Pi = points_S[:,i] + 1lines_S[:,i]
        y[i+N_poses] = norm(N_RB[:,i])^2-Ni'Ta
        A[i+N_poses,:] = (reshape(repeat(Ni'Ra,3,1)',9,1) .*[reshape(repeat(Pi[1:2],1,3)',6,1);1;1;1])'
    end
    A,y
end

function transform_data(POSES, points_S, lines_S, T_TF_S)
    N_poses = size(POSES,3)
    points  = zeros(4,N_poses)
    lines   = similar(lines_S)
    for i = 1:N_poses
        T_RB_S      = POSES[:,:,i] * T_TF_S
        points[:,i] = T_RB_S*[points_S[:,i]; 1]
        lines[:,i]  = T2R(T_RB_S)*lines_S[:,i]
    end
    points, lines
end

function findplanes(planes, points)
    N_poses = length(planes)
    N_planes  = maximum(planes)
    μ_RB    = zeros(3,N_poses)
    N_RB    = μ_RB
    normals = zeros(N_planes,3)
    for j = 1:N_planes
        ind          = planes .== j
        μᵢ           = mean(points[1:3,ind],dims=2)[:]
        μ_RB[:,ind]  = repeat(μᵢ,1,sum(ind))
        D,V          = eigen(cov(points[1:3,ind]')) # PCA
        N_RBi        = V[:,1] # Eigenvector of smalles eigval is the normal of the plane
        if μᵢ'N_RBi < 0 # Make sure all normals point in same direction
            N_RBi = -1*N_RBi
        end
        N_RBi        = N_RBi*(N_RBi'μᵢ)
        N_RB[:,ind]  = repeat(N_RBi,1,sum(ind)) # Repeat plane normal for all planes with same ind
        normals[j,:] = N_RBi
    end
    normals, N_RB
end

function plotLines(points,lines)
    P = size(points,2)
    for i = 1:P
        p = points[1:3,i]
        l = lines[1:3,i]
        plot3smart([p'+0.05*l'; p'-0.05*l'],c=:red)
    end
end

function plotPlanes(normals)
    P = size(normals,1)
    for i = 1:P
        n = normals[i,:]'
        xdir = n+[1, 0, 0]
        ydir = normalize(cross(n,xdir))
        xdir = normalize(cross(ydir,n))

        p[:,1] = n + xdir + ydir
        p[:,2] = n + xdir - ydir
        p[:,3] = n - xdir - ydir
        p[:,4] = n - xdir + ydir

        fill3(p[1,:],p[2,:],p[3,:])
    end
    plot!(alpha=0.5, zlims=(-0.3, 1))
end

function pointDiff(T_TF_S, T_RB_TF, points_S)
    N = size(points_S,2)
    @assert size(T_RB_TF,3) == N
    points_RB = zeros(4,N)
    for i = 1:N
        T_RB_S = T_RB_TF[:,:,i]*T_TF_S
        points_RB[:,i] = T_RB_S*[points_S[:,i];1]
    end
    points_RB = points_RB[1:3,:]
    S = eigvals(cov(points_RB'))
    d = minimum(S)
end
