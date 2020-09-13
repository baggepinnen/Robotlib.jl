using Robotlib, LinearAlgebra, Random, Statistics, Plots
using Test
import Robotlib: ad, adi, I4

function simulateCalibration1(N)
    n = 6
    q = 2π * rand(N, n)

    wn = zeros(3, n)
    pn = zeros(3, n)
    xin = zeros(6, n + 1)
    wn[:, 1] = [0, 0, 1]
    pn[:, 1] = [0, 0, 0]
    wn[:, 2] = [0, 1, 0]
    pn[:, 2] = [0, 0, 0.3]
    wn[:, 3] = [0, 1, 0]
    pn[:, 3] = [0, 0, 1]
    wn[:, 4] = [1, 0, 0]
    pn[:, 4] = [1, 0, 1]
    wn[:, 5] = [0, 1, 0]
    pn[:, 5] = [1, 0, 1]
    wn[:, 6] = [1, 0, 0]
    pn[:, 6] = [1, 0, 1]
    T0 = [
    0 -1 0 1.2
    0 0 -1 0
    1 0 0 1
    0 0 0 1
    ]
    for i = 1:n
        xin[:, i] = [-skew(pn[:, i]) * wn[:, i]; wn[:, i]]
    end
    xin[:, n+1] = twistcoords(logT(T0))
    xinmod = deepcopy(xin)
    for i = 1:n+1
        xinmod[:, i] += [0.001randn(3); 0.01π / 180 * randn(3)]
        #conformize(xinmod[:,i])
    end
    Ta = zeros(4, 4, N)
    for i = 1:N
        Ta[:, :, i] = fkinePOE(xin, q[i, :])
    end
    return q, xin, T0, xinmod, Ta
end

@testset "Robotlib tests" begin
    Random.seed!(1)

    @testset "calibNAXP" begin
        @info "Testing calibNAXP"
        include("testCalibNAXP.jl")

    end
    @testset "Force" begin
        @info "Testing Force"
        include("testForce.jl")

    end

    @testset "DH" begin
        @info "Testing DH"
        DH2400()
        DH7600()
        DHYuMi()
        dh = DHtest()
        DH2twistsPOE(dh)
        DH2twistsLPOE(dh)

    end

    @testset "Plotting" begin
        @info "Testing Plotting"
        q, xin, T0, xinmod, Ta = simulateCalibration1(10)
        trajplot(Ta)
        trajplot3(Ta)
    end


    @testset "Utils" begin
        @info "Testing Utils"
        R = toOrthoNormal(randn(3,3))
        @test Rangle(R,R) < 1e-7
        R2 = rpy2R(1*pi/180,0,0)*R
        @test Rangle(R,R2,true) <= 1.0000001

        q = Quaternion(R)
        @test Robotlib.Quaternions.rotationmatrix(q) ≈ R

        r,p,y = R2rpy(R, conv="xyz")
        @test rpy2R(r,p,y, "xyz") ≈ R

    end



    @testset "Calibration and kinematics" begin



        @testset "Inverse kinematics" begin
            @info "Testing Inverse kinematics"

            dh = DH7600();
            q0 = ones(6)
            xi = DH2twistsPOE(dh)
            T0 = fkinePOE(xi,q0)
            Tt = deepcopy(T0)
            Tt[1:3,4] += 0.1*[1,1,1]
            Tt[1:3,1:3] *= rpy2R(pi/180,0,0)

            @time q,err = ikinePOE(xi,Tt,q0,verbose=true)
            @test err[end] < sqrt(eps())

            # pp = semilogy(err)
            # display(pp)

        end

        @testset "Jacobian" begin
            @info "Testing Jacobian"

            dh = DH7600();
            q = [1,1,0,0,0,0];
            q = randn(6)
            AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, I4);
            Jn,J0,Tjac,Tijac,Tnjac = jacobian(q,dh, I4);

            xiL = DH2twistsLPOE(Tn0)
            xi = DH2twistsPOE(Tn0)
            Jsb, Jbb = jacobianPOE(q,xi)

            @test Jsb ≈ J0
            @test_broken Jbb ≈ Jn

            # display(round.(J0, digits=3))
            # display(round.(Jsb, digits=3))
        end


        N = 100

        # println("===== Testing calibPOE =====")
        # q, xin,T0, xinmod, Ta= Robotlib.Calibration.simulateCalibration_POE(N)
        # xic,et,er = Robotlib.Calibration.calibPOE(xinmod,Ta,q,maxiter=150, λ = 100.0)
        # println("Initial error: ", norm(xin[:,1:7]-xinmod[:,1:7]))
        # println("Final error: ", norm(xin[:,1:7]-xic[:,1:7]))
        # @test et[et .!= 0][end] < 1e-12
        # @test et[et .!= 0][end] < 1e-12

        println("===== Testing calibLPOE =====")
        q, xin, Tn0, Tn0mod, Ta =
            Robotlib.Calibration.simulateCalibration_LPOE(N)
        @time Tn0c, xic, et, er = Robotlib.Calibration.calibLPOE(
            xin,
            Tn0mod,
            Ta,
            q,
            maxiter = 10,
            λ = 100.0,
        )
        @test et[et.!=0][end] < 1e-12
        @test er[er.!=0][end] < 1e-12
        println("Initial error: ", sqrt(sum((Tn0 - Tn0mod) .^ 2)))
        println("Final error: ", sqrt(sum((Tn0 - Tn0c) .^ 2)))
        @test sqrt(sum((Tn0 - Tn0c) .^ 2)) < 0.8sqrt(sum((Tn0 - Tn0mod) .^ 2))

        println("===== Testing calibLPOEdual =====")
        q, xin, Tn0, Tn0mod, Ta =
            Robotlib.Calibration.simulateCalibration_LPOE_dual(N)
        @time Tn0c, xic, et, er = Robotlib.Calibration.calibLPOEdual(
            xin,
            Tn0mod,
            q,
            maxiter = 6,
            λ = 0.01,
        )
        println(
            "Error between Tn0 and Tn0mod: ",
            sqrt(
                sum(
                    (Tn0[1:3, 4, :] - Tn0mod[1:3, 4, :, 1]) .^ 2 +
                    (Tn0[1:3, 4, :] - Tn0mod[1:3, 4, :, 2]) .^ 2,
                ) / N,
            ),
        )
        println(
            "Error between Tn0 and Tn0c  : ",
            sqrt(
                sum(
                    (Tn0[1:3, 4, :] - Tn0c[1:3, 4, :, 1]) .^ 2 +
                    (Tn0[1:3, 4, :] - Tn0c[1:3, 4, :, 2]) .^ 2,
                ) / N,
            ),
        )
        global ei = 0.0
        global ec = 0.0
        for i = 1:N
            global ei, ec
            T1 = fkineLPOE(Tn0mod[:, :, :, 1], xin[:, :, 1], q[i, :, 1])
            T2 = fkineLPOE(Tn0mod[:, :, :, 2], xin[:, :, 2], q[i, :, 2])
            ei += norm(twistcoords(log(Ta[:, :, i] * trinv(T1))))
            ei += norm(twistcoords(log(Ta[:, :, i] * trinv(T2))))
            T1 = fkineLPOE(Tn0c[:, :, :, 1], xin[:, :, 1], q[i, :, 1])
            T2 = fkineLPOE(Tn0c[:, :, :, 2], xin[:, :, 2], q[i, :, 2])
            ec += norm(twistcoords(log(Ta[:, :, i] * trinv(T1))))
            ec += norm(twistcoords(log(Ta[:, :, i] * trinv(T2))))
        end
        println(
            "Initial error: ",
            round(ei / N, digits = 5),
            " Calibrated error: ",
            round(ec / N, digits = 5),
        )
        @test ec < ei

    end


    @testset "Frames" begin
        @info "Testing Frames"
        include("testFrames.jl")
    end
end
