using Robotlib
using Base.Test
import Robotlib: ad, adi


function simulateCalibration1(N)
    n       = 6
    q       = 2π*rand(N,n)

    wn      = zeros(3,n)
    pn      = zeros(3,n)
    xin     = zeros(6,n+1)
    wn[:,1] = [0,0,1]
    pn[:,1] = [0,0,0]
    wn[:,2] = [0,1,0]
    pn[:,2] = [0,0,0.3]
    wn[:,3] = [0,1,0]
    pn[:,3] = [0,0,1]
    wn[:,4] = [1,0,0]
    pn[:,4] = [1,0,1]
    wn[:,5] = [0,1,0]
    pn[:,5] = [1,0,1]
    wn[:,6] = [1,0,0]
    pn[:,6] = [1,0,1]
    T0 = [  0 -1 0  1.2;
    0 0 -1 0;
    1 0  0  1;
    0 0  0  1]
    for i = 1:n
        xin[:,i] = [-skew(pn[:,i])*wn[:,i]; wn[:,i]]
    end
    xin[:,n+1] = twistcoords(logT(T0))
    xinmod   = deepcopy(xin)
    for i = 1:n+1
        xinmod[:,i]  += [0.001randn(3); 0.01π/180*randn(3)]
        #conformize(xinmod[:,i])
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkinePOE(xin,q[i,:]')
    end
    return q, xin, T0, xinmod, Ta
end

function simulateCalibration2(N)
    n       = 6
    q       = 2π*rand(N,n)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, eye(4));
    xin     = DH2twistsPOE(dh)
    xinmod   = deepcopy(xin)
    for i = 1:n+1
        xinmod[:,i]  += [0.001randn(3); 0.01π/180*randn(3)]
        #conformize(xinmod[:,i])
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkinePOE(xin,q[i,:]')
    end
    return q, xin, T0, xinmod, Ta
end

function simulateCalibrationLPOE(N)
    srand(1)
    n       = 6
    q       = 2π*rand(N,n)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, eye(4));
    xin = DH2twistsLPOE(Tn0)
    Tn0mod = deepcopy(Tn0)
    for i = 1:n+1
        Tn0mod[1:3,4,i] += 0.1randn(3)
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkineLPOE(Tn0,xin,q[i,:]')
        AAA,BBB,T = jacobian(q[i,:]',dh, eye(4));
        @assert T ≈ Ta[:,:,i]
    end
    return q, xin, Tn0, Tn0mod, Ta
end


function simulateCalibration4(N)
    #srand(1)
    n       = 6
    qt      = 2π*rand(N,n)
    q1      = qt + 0.01*π/180*randn(N,n)
    q2      = qt + 0.01*π/180*randn(N,n)
    q       = cat(3,q1,q2)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, eye(4));
    xin = DH2twistsLPOE(Tn0)
    xin = cat(3,xin,xin)
    Tn0mod = cat(4,deepcopy(Tn0),deepcopy(Tn0))
    for i = 1:n+1
        Tn0mod[1:3,4,i,1] += 0.1randn(3)
        Tn0mod[1:3,4,i,2] += 0.1randn(3)
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkineLPOE(Tn0,xin,qt[i,:]')
    end

    return q, xin, Tn0, Tn0mod, Ta
end

function run_tests()
    # q, xin,T0, xinmod, Ta= simulateCalibration2(100)
    # xic,et,er = Robotlib.Calibration.calibPOE(xinmod,Ta,q,maxiter=150, λ = 10000.0)
    # display(norm(xin[:,1:7]-xinmod[:,1:7]))
    # display(norm(xin[:,1:7]-xic[:,1:7]))
    N = 100
    println("===== Testing calibLPOE =====")
    q, xin,Tn0,Tn0mod, Ta = simulateCalibrationLPOE(N)
    @time Tn0c,xic,et,er = Robotlib.Calibration.calibLPOE(xin,Tn0mod,Ta,q,maxiter=10, λ = 100.0)
    @test et[et .!= 0][end] < 1e-12
    @test er[er .!= 0][end] < 1e-12
    println("Initial error: ", sqrt(sum((Tn0-Tn0mod).^2)))
    println("Final error: ", sqrt(sum((Tn0-Tn0c).^2)))
    @test sqrt(sum((Tn0-Tn0c).^2)) < 0.8sqrt(sum((Tn0-Tn0mod).^2))

    println("===== Testing calibLPOEdual =====")
    q, xin,Tn0,Tn0mod, Ta = simulateCalibration4(N)
    @time Tn0c,xic,et,er = Robotlib.Calibration.calibLPOEdual(xin,Tn0mod,q,maxiter=6, λ = 0.01)
    println("Error between Tn0 and Tn0mod: ",
        sqrt(sum((Tn0[1:3,4,:]-Tn0mod[1:3,4,:,1]).^2+(Tn0[1:3,4,:]-Tn0mod[1:3,4,:,2]).^2)/N))
    println("Error between Tn0 and Tn0c  : ",
        sqrt(sum((Tn0[1:3,4,:]-Tn0c[1:3,4,:,1]).^2+(Tn0[1:3,4,:]-Tn0c[1:3,4,:,2]).^2)/N))
    ei = 0.0
    ec = 0.0
    for i = 1:N
        T1 = fkineLPOE(Tn0mod,xin,q[i,:,1]')
        T2 = fkineLPOE(Tn0mod,xin,q[i,:,2]')
        ei += norm(twistcoords(logm(Ta[:,:,i]*trinv(T1))))
        ei += norm(twistcoords(logm(Ta[:,:,i]*trinv(T2))))
        T1 = fkineLPOE(Tn0c,xin,q[i,:,1]')
        T2 = fkineLPOE(Tn0c,xin,q[i,:,2]')
        ec += norm(twistcoords(logm(Ta[:,:,i]*trinv(T1))))
        ec += norm(twistcoords(logm(Ta[:,:,i]*trinv(T2))))
    end
    println("Initial error: ",round(ei/N,5), " Calibrated error: ", round(ec/N,5))
    @test ec < ei
end

run_tests()

println("===== Testing frames =====")
include("testFrames.jl")
