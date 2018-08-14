using Robotlib
import Robotlib: Rt2T, T2R, T2t
include(joinpath(dirname(@__FILE__),"..","src","calibNAXP.jl"))

mutable struct Plane
    N
    N1
    d
    Pn
    Qn
end

function generateRandomPlane(no)
    if no == 1
        N_RB = [1, 0, 0]
    elseif no == 2
        N_RB = [0, 1, 0]
    elseif no == 3
        N_RB = [0, 0, 1]
    else
        N_RB = normalize(randn(3)) # Normal of calibration plane, given in RB
    end
    d_RB = 1 #abs(0.001*randn);
    Pn = N_RB*N_RB' # Projection matrix onto normal
    Qn = eye(3)-Pn # Projection matrix onto plane

    Plane(N_RB, [N_RB;1], d_RB, Pn, Qn)
end

function generateRandomPoses(N_poses, sigma_n,T_TF_S, plane)
    LN_S    = [0, 0, 1] # Normal of laser plane in sensor frame
    N_RB    = plane.N
    d_RB    = plane.d
    T_RB_TF = Array(Float64,4,4,N_poses)
    i       = 1
    lines_RB= zeros(3,N_poses)
    lines_S = zeros(3,N_poses)
    points_RB = zeros(3,N_poses)
    points_S = zeros(4,N_poses)
    while i <= N_poses
        T_RB_TF[:,:,i]      = Rt2T(rpy2R(1.5*randn(), 1.5*randn(), 1.5*randn()),0.6*randn(3) + 2*N_RB)
        T_RB_S              = T_RB_TF[:,:,i]*T_TF_S
        LN_RB               = T2R(T_RB_S)*LN_S
        intersection_line_RB= normalize(cross(N_RB,LN_RB))
        intersection_line_S = inv(T2R(T_RB_S))*intersection_line_RB
        lines_RB[:,i]       = intersection_line_RB
        lines_S[:,i]        = intersection_line_S
        laser_0_RB          = T_RB_S[1:3,4]
        laser_y_RB          = T_RB_S[1:3,2]
        lambda              = (d_RB-vecdot(N_RB,laser_0_RB))/(vecdot(N_RB,laser_y_RB))
        points_RB[:,i]      = laser_0_RB + lambda*laser_y_RB
        points_S[:,i]       = inv(T_RB_S)*[points_RB[:,i]; 1]
        if abs(lambda) > 2
            continue
        end
        i += 1
    end
    points_S = points_S[1:3,:]
    points_S = points_S + sigma_n*randn(size(points_S))

    # Verify poses
    C = cov(points_RB');    D,V = eig(C)
    nhat_points = V[:,1]
    C = cov(lines_RB');    D,V = eig(C)
    nhat_lines = V[:,1]
    if 1-abs(vecdot(nhat_points,nhat_lines)) > 0.001 || D[1,1] > 1e-5
        @warn("Something is wrong")
        @show 1-abs(vecdot(nhat_points,nhat_lines)), D[1,1]
    end

    points = zeros(4,N_poses)
    for i = 1:N_poses
        T_RB_S = T_RB_TF[:,:,i] * T_TF_S
        #         points(:,end+1) = T_RB_S*[(points_S[:,i]);1];
        points[:,i] = T_RB_S*[(points_S[:,i]);1]
        #         points(:,end+1) = inv(T_RB_S)*[(points_S[:,i] - 0.01*lines_S[:,i]);1];
    end
    if any(abs(points[1:3,:] - points_RB) .> 1e-2)
        @warn("Something is wrong 2")
    end
    points_S,lines_S, T_RB_TF
end


Random.seed!(1)
MC = 100
N_planes = 3
N_points = 2
N_calibs = 30
N_poses = 10
sigma_n = 5e-4
SSEStart = zeros(MC)
normStart = zeros(MC)
distStart = zeros(MC)
rotStart = zeros(MC)
SSEMC = zeros(MC,N_calibs)
normMC = zeros(MC,N_calibs)

for mc = 1:MC
    # Simulate measurement data
    T_TF_S = Rt2T(rpy2R(randn(3)),0.5randn(3)) #TCP to sensor
    points_S = Array(Float64,3,0);    lines_S = Array(Float64,3,0);    T_RB_TF = Array(Float64,4,4,0);    plane = Plane[]; planes = Array(Int,0)
    SSE = zeros(N_planes)
    for j = 1:N_planes
        push!(plane, generateRandomPlane(j))
        points_St,lines_St, T_RB_TFt = generateRandomPoses(N_poses, sigma_n,T_TF_S, plane[j])
        points_S = [points_S points_St]
        lines_S = [lines_S lines_St]
        T_RB_TF = cat(3,T_RB_TF, T_RB_TFt)
        planes = [planes; j*ones(Int,N_poses)]
        SSE[j] = pointDiff(T_TF_S,T_RB_TFt,points_St)[1]
    end
    @edb any(abs(SSE .> 1e-3)) && warn("Points does not seem to lie on a plane")
    dist = [norm(T_RB_TF[1:4,4,i] - T_TF_S*[points_S[:,i];1]) for i in 1:size(T_RB_TF,3)]
    #     hist(dist,15)

    ## Calibration
    # display('----------------------------------------')
    # Generate nominal transformation
    N_poses = size(points_S,2)
    T_TF_S_real = T_TF_S
    if true
        T_TF_S0 = Rt2T(rpy2R(60π/180*(rand(3)-0.5)),0.4*(rand(3)-0.5))*T_TF_S_real
        #         T_TF_S0 = rt2tr(rpy2R(2*(rand(1,3)-0.5),'deg'),0.01*(rand(3,1)-0.5))*T_TF_S_real;
    else
        T_TF_S0 = T_TF_S_real
    end

    for j = 1:N_planes
        ind = planes .== j
        ind = findall(ind)
        SSE[j] = sqrt(pointDiff(T_TF_S0,T_RB_TF[:,:,ind],points_S[1:3,ind])[1])
    end

    SSEStart[mc] = mean(SSE)
    normStart[mc] = vecnorm(T_TF_S0-T_TF_S_real)
    distStart[mc] = norm(T_TF_S0[1:3,4]-T_TF_S_real[1:3,4])
    rotStart[mc] = norm(180/π*R2rpy(T_TF_S0\T_TF_S_real))


    T_TF_S, RMScalibs, norms = calibNAXP(points_S, lines_S, T_RB_TF, T_TF_S0, planes,  N_calibs)#, T_TF_S_real)
    SSEMC[mc,:] = RMScalibs
    normMC[mc,:] = norms

    if RMScalibs[end] > 1e-2
        println("Bad result")
    end

    distEnd[mc] = norm(T_TF_S[1:3,4]-T_TF_S_real[1:3,4])
    rotEnd[mc] = norm(180/π*R2rpy(T_TF_S\T_TF_S_real))
end
# normStart
# normMC


plot(0:N_calibs,[normStart', normMC]',yscale=:log10,c=:black, xlabel="Number of iterations")
plot!([0 N_calibs], [sigma_n sigma_n],l=:dash, c=:red)

plot(0:N_calibs,[SSEStart', SSEMC]',yscale=:log10,c=:black, xlabel="Number of iterations",title="RMS distance from points to plane [m]")
plot!([0 N_calibs], [sigma_n sigma_n],l=:dash,c=:red)

boxplot(["Before" "After"],([distStart', distEnd']), title="Distance error [m]")
# hold on, plot([0 3], [sigma_n sigma_n],'--r')

boxplot(["Before" "After"],([rotStart', rotEnd']),title="Rotation error [degree]")
