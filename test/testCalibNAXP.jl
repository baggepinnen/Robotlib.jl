using Plots
using Statistics, LinearAlgebra, Random
Random.seed!(1)
using Robotlib, StaticArrays, Test
import Robotlib: Rt2T, T2R, T2t, I4
# include(joinpath(dirname(@__FILE__),"..","src","calibNAXP.jl"))

const TM = SMatrix{4,4,Float64,16}

mutable struct Plane{T1,T2,T3,T4,T5}
    N::T1
    N1::T2
    d::T3
    Pn::T4
    Qn::T5
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
    d_RB = 1          # abs(0.001*randn);
    Pn = N_RB * N_RB' # Projection matrix onto normal
    Qn = I - Pn       # Projection matrix onto plane

    Plane(N_RB, [N_RB; 1], d_RB, Pn, Qn)
end

function generateRandomPoses(N_poses, σₙ, T_TF_S, plane)
    LN_S      = [0, 0, 1] # Normal of laser plane in sensor frame
    N_RB      = plane.N
    d_RB      = plane.d
    T_RB_TF   = Array{Float64}(undef, 4, 4, N_poses)
    lines_RB  = zeros(3, N_poses)
    lines_S   = zeros(3, N_poses)
    points_RB = zeros(3, N_poses)
    points_S  = zeros(4, N_poses)
    i = 1
    while i <= N_poses
        T_RB_TF[:, :, i] = Rt2T(
            rpy2R(1.5 * randn(), 1.5 * randn(), 1.5 * randn()),
            0.6 * randn(3) + 2 * N_RB,
        )
        T_RB_S               = T_RB_TF[:, :, i] * T_TF_S |> TM
        LN_RB                = T2R(T_RB_S) * LN_S
        intersection_line_RB = normalize(cross(N_RB, LN_RB))
        intersection_line_S  = inv(T2R(T_RB_S)) * intersection_line_RB
        lines_RB[:, i]       = intersection_line_RB
        lines_S[:, i]        = intersection_line_S
        laser_orig_RB        = T_RB_S[1:3, 4]
        laser_y_RB           = T_RB_S[1:3, 2]
        λ                    = (d_RB - N_RB'laser_orig_RB) / (N_RB'laser_y_RB)
        # λ is the distance from laser origin to plane, should be positive and small to be realistic, also, the angle with the plane normal should not be too large. In paper, the following (slightly unrealistic) condition was used
        # if abs(λ) > 2
        if abs(λ) > 0.2 || λ < 0 || (180 / pi .* acos(intersection_line_S[2])) > 120
            continue
        end
        points_RB[:, i] = laser_orig_RB + λ * laser_y_RB
        points_S[:, i] = inv(T_RB_S) * [points_RB[:, i]; 1]
        i += 1
    end
    points_S = points_S[1:3, :]
    points_S .+= σₙ * randn(size(points_S))

    # Verify poses
    ei          = eigen(cov(points_RB'))
    nhat_points = ei.vectors[:, 1]
    V           = eigen(cov(lines_RB')).vectors
    nhat_lines  = V[:, 1]
    if 1 - abs(nhat_points'nhat_lines) > 0.001 || ei.values[1, 1] > 1e-5
        @warn("Something is wrong")
        @show 1 - abs(nhat_points'nhat_lines), ei.values[1, 1]
    end

    points = zeros(4, N_poses)
    for i = 1:N_poses
        T_RB_S = T_RB_TF[:, :, i] * T_TF_S |> TM
        #         points(:,end+1) = T_RB_S*[(points_S[:,i]);1];
        points[:, i] = T_RB_S * [points_S[:, i]; 1]
        #         points(:,end+1) = inv(T_RB_S)*[(points_S[:,i] - 0.01*lines_S[:,i]);1];
    end
    if any(abs.(points[1:3, :] - points_RB) .> 1e-2)
        @warn("Something is wrong 2")
    end
    points_S, lines_S, T_RB_TF
end



MC = 10   # In paper 100
N_planes = 3    # In paper 3
iters = 80   # In paper 30
N_poses = 20   # In paper 10
σₙ = 5e-4 # In paper: 5e-4
function run_calib(verbose = false)
    SSEStart  = zeros(MC)
    normStart = zeros(MC)
    distStart = zeros(MC)
    rotStart  = zeros(MC)
    distEnd   = zeros(MC)
    rotEnd    = zeros(MC)
    SSEMC     = zeros(MC, iters)
    normMC    = zeros(MC, iters)
    distEnd2  = zeros(MC)
    rotEnd2   = zeros(MC)
    SSEMC2    = zeros(MC, iters)
    normMC2   = zeros(MC, iters)
    mc = 1
    for mc = 1:MC
        @info "MC iterations $(mc)/$(MC)"
        # Simulate measurement data
        T_TF_S = Rt2T(rpy2R(randn(3)), 0.5 * randn(3)) |> TM #TCP to sensor
        points_Sv = Matrix{Float64}[]
        lines_Sv = Matrix{Float64}[]
        T_RB_TFv = Array{Float64,3}[]
        plane = Plane[]
        SSE = zeros(N_planes)
        verbose && println("Genereating planes and poses")
        for j = 1:N_planes
            push!(plane, generateRandomPlane(j))
            points_St, lines_St, T_RB_TFt =
                generateRandomPoses(N_poses, σₙ, T_TF_S, plane[j])
            push!(points_Sv, points_St)
            push!(lines_Sv, lines_St)
            push!(T_RB_TFv, T_RB_TFt)
            SSE[j] = Robotlib.pointDiff(T_TF_S, T_RB_TFt, points_St)[1]
        end
        planes = repeat((1:N_planes)', N_poses)[:]
        points_S = cat(points_Sv..., dims = 2)
        lines_S = cat(lines_Sv..., dims = 2)
        T_RB_TF = cat(T_RB_TFv..., dims = 3)

        any(abs.(SSE .> 1e-3)) && @error("Points does not seem to lie on a plane")
        dist = [
            norm(T_RB_TF[1:4, 4, i] - T_TF_S * [points_S[:, i]; 1])
            for i = 1:size(T_RB_TF, 3)
        ]
        #     hist(dist,15)

        ## Calibration
        # display('----------------------------------------')
        # Generate nominal transformation
        T_TF_S_real = copy(T_TF_S)
        if true
            T_TF_S0 =
                Rt2T(rpy2R(60π / 180 * (rand(3) .- 0.5)), 0.4 * (rand(3) .- 0.5)) *
                T_TF_S_real
            # T_TF_S0 = Rt2T(rpy2R(10π/180*(rand(3).-0.5)),0.1*(rand(3).-0.5))*T_TF_S_real
            #         T_TF_S0 = rt2tr(rpy2R(2*(rand(1,3).-0.5),'deg'),0.01*(rand(3,1).-0.5))*T_TF_S_real;
        else
            T_TF_S0 = copy(T_TF_S_real)
        end

        for j = 1:N_planes
            ind = findall(planes .== j)
            SSE[j] = sqrt(Robotlib.pointDiff(T_TF_S0, T_RB_TF[:, :, ind], points_S[1:3, ind])[1])
        end

        SSEStart[mc] = mean(SSE)
        normStart[mc] = norm(vec(T_TF_S0 - T_TF_S_real))
        distStart[mc] = norm(T_TF_S0[1:3, 4] - T_TF_S_real[1:3, 4])
        rotStart[mc] = norm(180 / π * R2rpy(T_TF_S0 \ T_TF_S_real))
        # display(T_TF_S_real)
        T_TF_S, RMScalibs, ALLcalibs, norms = calibNAXP(
            points_S,
            lines_S,
            T_RB_TF,
            T_TF_S0,
            planes,
            iters,
            trueT_TF_S = T_TF_S_real,
            variance = (1, 1),
            doplot = true,
        )
        T_TF_S2, RMScalibs2, ALLcalibs2, norms2 = T_TF_S, RMScalibs, ALLcalibs, norms#Robotlib.Calibration.calibNAXP_bootstrap(points_S, lines_S, T_RB_TF, T_TF_S0, planes,  iters, trueT_TF_S=T_TF_S_real, nbootstrap=500)

        SSEMC[mc, :] = RMScalibs
        normMC[mc, :] = norms
        if RMScalibs[end] > 1e-2
            println("Bad result")
        end
        distEnd[mc] = norm(T_TF_S[1:3, 4] - T_TF_S_real[1:3, 4])
        rotEnd[mc] = norm(180 / π * R2rpy(T_TF_S \ T_TF_S_real))

        SSEMC2[mc, :] = RMScalibs2
        normMC2[mc, :] = norms2
        if RMScalibs2[end] > 1e-2
            println("Bad result")
        end
        distEnd2[mc] = norm(T_TF_S2[1:3, 4] - T_TF_S_real[1:3, 4])
        rotEnd2[mc] = norm(180 / π * R2rpy(T_TF_S2 \ T_TF_S_real))
    end
    SSEStart,
    normStart,
    distStart,
    rotStart,
    distEnd,
    rotEnd,
    SSEMC,
    normMC,
    distEnd2,
    rotEnd2,
    SSEMC2,
    normMC2
end

# @testset "CalibNAXP" begin

@test isapprox(Robotlib.pointDiff(I4, cat(fill(I4, 5)..., dims = 3), zeros(3, 5))[], 0, atol = 1e-10)

SSEStart, normStart, distStart, rotStart, distEnd, rotEnd, SSEMC, normMC, distEnd2, rotEnd2, SSEMC2, normMC2 = run_calib()

# iters = size(SSEMC,2)
# using StatPlots
# gr(legend=false)
# plot(0:iters,copy([normStart normMC]'),yscale=:log10,c=:black, xlabel="Number of iterations", layout=2, subplot=1)
# # plot!(0:iters,copy([normStart normMC2]'),yscale=:log10,c=:green, subplot=1)
# hline!([σₙ, σₙ],l=:dash, c=:red, subplot=1)
#
# plot!(0:iters,[SSEStart SSEMC]',yscale=:log10,c=:black, xlabel="Number of iterations",title="RMS distance from points to plane [m]", subplot=2)
# # plot!(0:iters,[SSEStart SSEMC2]',c=:green, subplot=2)
# hline!([σₙ σₙ],l=:dash,c=:red, subplot=2)
#
#
# boxplot(["Before" "After" "After WTLS"],([distStart distEnd distEnd2]), title="Distance error [m]", yscale=:log10, layout=2, subplot=1)
# hline!([σₙ σₙ],c=:red)
# boxplot!(["Before" "After" "After WTLS"],([rotStart rotEnd rotEnd2]),title="Rotation error [degree]", yscale=:log10, subplot=2)

# end
