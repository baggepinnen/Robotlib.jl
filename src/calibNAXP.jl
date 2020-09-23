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
function calibNAXP(points_S, lines_S, POSES, T_TF_S, planes::AbstractVector{Int},  iters::Integer = 50; doplot=false, trueT_TF_S=nothing, variance=nothing)
    N_poses   = size(POSES,3)
    N_planes  = maximum(planes)
    RMScalibs = zeros(iters)
    ALLcalibs = zeros(iters, N_planes)
    norms     = trueT_TF_S == nothing ? nothing : zeros(iters)
    local points, lines, N_RB
    for c = 1:iters
        # Convert points to RB cordinate system
        points, lines = transform_data(POSES, points_S, lines_S, T_TF_S)
        N_RB = findplanes(planes, points)

        # Estimation
        if variance == nothing
            A,y = get_matrices(N_RB, POSES, points_S, lines_S)
            w   = A\y
        else
            A,y,Σy = get_matrices(N_RB, POSES, points_S, lines_S, variance)
            w   = wls(A,y,Σy)
        end
        # w   = tls(A,y)
        R,t = w2Rt(w)

         # Re-estimation of translation
        w₂     = R[1:3,1:2][:]
        y₂     = y - A[:,1:6]*w₂
        t₂     = A[:,7:9]\y₂
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
        try
            plotPlanes(unique(N_RB,dims=2))
            plotLines(points,lines)
        catch ex
            if ex isa UndefVarError
                @error "Plotting failed. Plotting is only available if you have manually run `using Plots`"
            else
                rethrow(ex)
            end
        end
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
    R   = orthonormal(Rnn)
    t   = H[1:3,3]
    R,t
end

function get_matrices(N_RB, POSES, points_S, lines_S,variance=nothing)
    N_poses = size(POSES,3)
    A = zeros(2*N_poses,9)
    y = zeros(2*N_poses)
    𝜸 = 1
    if variance != nothing
        α = mapslices(lines_S, dims=1) do l
            acos(l[2]) # = acos(eʸ'l)
        end[:]
        Σ(α) = [sin(α)^2 -sin(α)*cos(α); -sin(α)*cos(α) cos(α)^2]
        Σp,σα = variance
        Σp2 = [Σp .+ 𝜸^2*σα^2*abs.(Σ(α)) for α in α]
        Σy = spzeros(2*N_poses,2*N_poses)
    end
    for i = 1:N_poses
        Ni     = N_RB[:,i]
        Ra     = POSES[1:3,1:3,i]
        Ta     = POSES[1:3,4,i]
        Pi     = points_S[:,i]
        y[i]   = norm(Ni)^2-Ni'Ta
        A[i,:] = (reshape(repeat(Ni'Ra,3,1)',9,1).*[reshape(repeat(Pi[1:2],1,3)',6,1);1;1;1])' #TODO: rewrite
        # Sample one additional point along line
        Pi = points_S[:,i] + 𝜸*lines_S[:,i]
        y[i+N_poses] = norm(N_RB[:,i])^2-Ni'Ta
        A[i+N_poses,:] = (reshape(repeat(Ni'Ra,3,1)',9,1) .*[reshape(repeat(Pi[1:2],1,3)',6,1);1;1;1])'
        if variance != nothing
            Σy[i,i] = Σp
            Σy[i+N_poses, i+N_poses] = tr(Σp2[i])
            Σy[i+N_poses, i] = Σy[i, i+N_poses] = Σp
        end
    end
    variance == nothing ? (A,y) : (A,y,Symmetric(Σy))
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

function findplanes(planes, points, rep=true)
    N_poses  = length(planes)
    N_planes = maximum(planes)
    N_out    = rep ? N_poses : N_planes
    μ_RB     = zeros(3,N_out)
    N_RB     = μ_RB
    for j = 1:N_planes
        pointind           = planes .== j
        planeind           = rep ? pointind : j

        μᵢ                 = mean(points[1:3,pointind],dims = 2)[:]
        μ_RB[:,planeind]  .= μᵢ
        D,V                = eigen(cov(points[1:3,pointind]')) # PCA
        N_RBi              = V[:,1] # Eigenvector of smalles eigval is the normal of the plane
        if μᵢ'N_RBi < 0 # Make sure all normals point in same direction
            N_RBi = -1*N_RBi
        end
        N_RBi             = N_RBi*(N_RBi'μᵢ)
        N_RB[:,planeind] .= N_RBi # Repeat plane normal for all planes with same ind
    end
    N_RB
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


## The code below is used to verify the statements made in the appendix to chapter 13 in the PhD thesis cited in the README.

# function bootstrap(fun::Function, estimate, N, k=100)
#     k < 1 && return (fun(1:N), nothing)
#     results = map(1:k) do _
#         sample  = rand(1:N, N)
#         result  = fun(sample)
#     end
#     fun(1:N), estimate(results)
# end
#
#
# function columnwise_cov(data)
#     numcols = size(data[1], 2)
#     map(1:numcols) do c
#         columndata = hcat([d[:,c] for d in data]...)'
#         cov(columndata)
#     end
# end
#
# """
# This function is the same as `calibNAXP` except it estimates the covariance of the estimated normals using bootstrapping, and solves the optimization problem using the weighted total least-squares algorithm instead of standard least-square. It does not really provide any additional benefints over the standard algorithm, but is considerable slower.
# """
# function calibNAXP_bootstrap(points_S, lines_S, POSES, T_TF_S, planes::AbstractVector{Int},  iters::Integer = 50; trueT_TF_S=nothing, nbootstrap=500)
#     N_poses   = size(POSES,3)
#     N_planes  = maximum(planes)
#     RMScalibs = zeros(iters)
#     ALLcalibs = zeros(iters, N_planes)
#     norms     = trueT_TF_S == nothing ? nothing : zeros(iters)
#     for c = 1:iters
#         # Convert points to RB cordinate system
#         points, lines = transform_data(POSES, points_S, lines_S, T_TF_S)
#         N_RB, Σn = Robotlib.Calibration.bootstrap(columnwise_cov, N_poses, nbootstrap) do sample
#             findplanes(planes[sample], points[:,sample], false)
#         end
#         N_RB = N_RB[:,planes]
#         Σn = Σn[planes]
#         rowQ = map(1:N_poses) do i
#             Q,T,p,n = Σn[i], POSES[:,:,i], points_S[:,i], N_RB[:,i]
#             R = T[1:3,1:3]
#             t = T[1:3,4]
#             Qaa = R'Q*R
#             σ = [p[1], p[2], 1]*[p[1], p[2], 1]'
#             QAA = kron(ones(3,3), Qaa) .* kron(σ, ones(3,3))
#             Qyy = t'Q*t # TODO: + cov of norm(n)^2 Q
#             Qay = kron(ones(3), R'Q*(n-t)) .* kron([p[1], p[2], 1], ones(3))
#
#             [QAA Qay; [Qay' Qyy]]
#         end
#         Qaa,Qay,Qyy = rowcovariance(rowQ)
#         # Estimation
#         A,y = get_matrices(N_RB, POSES, points_S, lines_S)
#         w   = wtls(A,y,kron([1 0.5;0.5 1],Qaa),kron([1 0.5;0.5 1],Qay),kron([1 0.5;0.5 1],Qyy), iters=15)
#         # w   = tls(A,y)
#         R,t = w2Rt(w)
#
#          # Re-estimation of translation
#         w₂     = R[1:3,1:2][:]
#         y₂     = y - A[:,1:6]*w₂
#         t₂     = c <= 5 ? A[:,7:9]\y₂ : tls(A[:,7:9],y₂)
#         T_TF_S = Rt2T(R,t₂)
#
#         # Calculate iteration errors
#         RMSi = map(1:N_planes) do j
#             ind     = planes .== j
#             ind     = findall(ind)
#             sqrt(pointDiff(T_TF_S,POSES[:,:,ind],points_S[1:3,ind])[1])
#         end
#         RMScalibs[c]   = mean(RMSi)
#         ALLcalibs[c,:] = RMSi
#         trueT_TF_S == nothing || (norms[c] = norm(T_TF_S-trueT_TF_S))
#         #  any(abs(RMSi) > 1e-1) && warn("Points does not seem to lie on a plane")
#     end
#     T_TF_S, RMScalibs,ALLcalibs, norms
# end
