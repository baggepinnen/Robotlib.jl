import Robotlib: ad, adi
using Random
"""
    Tn0, xi, et, er = calibLPOE(xin,Tn0in,Ta,q; maxiter=10, λ=1.0)

Performs knematic calibration using the local POE formulation of kinematics.
- `Tn0` Nominal forward kinematics
- `xi` Twists
- `et` Translational errors as function of iteration of algorithm
- `er` Rotational errors as function of iteration of algorithm
- `λ` Regularization parameter

See `runtests.jl` for simulated usage. Some functions to generate simulation data for calibration purposes are provided in $(@__FILE__)
"""
function calibLPOE(xin,Tn0in,Ta,q;maxiter=10, λ=1.0)
    xi      = deepcopy(xin)
    Tn0     = deepcopy(Tn0in)
    n       = size(xi,2)-1
    N       = size(q,1)
    y       = zeros(6N)
    A       = zeros(6N,6*(n+1))
    xini    = deepcopy(xin)
    Tn0cand = deepcopy(Tn0)
    v       = zeros(3) # only revolute joints
    et      = zeros(maxiter+1)
    er      = zeros(maxiter+1)
    et[1],er[1]   = evalError(xi,Tn0,Ta,q)
    println("Error: ",et[1], " λ: ", λ)
    @assert size(Ta,3) == N

    for iter = 1:maxiter # do a few iterations of the calibration
        for i = 1:N
            Tfull = fkineLPOE(Tn0,xi,q[i,:])
            # populate the matrices of the linear estimation problem
            y[(i-1)*6+1:6i] = twistcoords(logT(Ta[:,:,i]*trinv(Tfull)))
            A[(i-1)*6+1:6*i,1:6] = ad(Tn0[:,:,1])
            adi = ad(Tn0[:,:,1]*expξ(xi[:,1],q[i,1]))
            for j = 2:n
                A[(i-1)*6+1:6*i,(j-1)*6+1:6*j] = adi*ad(Tn0[:,:,j])
                adi = adi*ad(Tn0[:,:,j]*expξ(xi[:,j],q[i,j]))
            end
            A[(i-1)*6+1:6*i,n*6+1:6*(n+1)] = adi*ad(Tn0[:,:,n+1])
        end
        x = (A'A + λ*I)\A'y
        # Update the candidate nominal parameters
        for j = 1:n+1
            ξ = x[(j-1)*6+1:6j]
            Tn0cand[:,:,j] = Tn0[:,:,j]*expξ2(ξ,1)
        end

        et[iter+1],er[iter+1] = evalError(xi,Tn0cand,Ta,q)
        println("Error: ",round(et[iter+1], digits=5), " λ: ", λ, " Norm dx: ", round(norm(x), digits=5))

        if norm(x) < 1e-16
            return Tn0,xi, et, er
        end
        if et[iter+1] > et[iter]
            et[iter+1] = et[iter]
            er[iter+1] = er[iter]
            λ *= 10
        else
            λ/=10
            xi = deepcopy(xini)
            Tn0 = deepcopy(Tn0cand)
        end
    end
    return Tn0, xi, et, er
end


"""
`calibLPOEdual(xi,Tn0,q;maxiter=10, λ=1.0)`
Performs dual arm calibration given joint twists xi and nominal transformations Tn0
Returns calibrated nominal transformations.
This function does not converge to the correct kinematic parameters since the problem is underspecified. It only converges to kinematic parameters that make the two arms agree.
"""
function calibLPOEdual(xin,Tn0in,q;maxiter=10, λ=1.0)

    xi      = deepcopy(xin)
    Tn0     = deepcopy(Tn0in)
    n       = size(xi,2)-1
    N       = size(q,1)
    y       = zeros(6N)
    A       = zeros(6N,6*(n+1),2)
    xini    = deepcopy(xin)
    Tn0cand = deepcopy(Tn0)
    et      = zeros(maxiter+1)
    er      = zeros(maxiter+1)
    et[1],er[1]   = evalErrorDual(xi,Tn0,q)
    println("Error: ",et[1], " λ: ", λ)

    for iter = 1:maxiter # do a few iterations of the calibration
        for i = 1:N
            T1 = fkineLPOE(Tn0[:,:,:,1],xi[:,:,1],q[i,:,1])
            T2 = fkineLPOE(Tn0[:,:,:,2],xi[:,:,2],q[i,:,2]) # Use the other arm as measurement
            y[(i-1)*6+1:6i] = twistcoords(logT(T2*trinv(T1)))
            for arm = 1:2
                A[(i-1)*6+1:6*i,1:6,arm] = ad(Tn0[:,:,1,arm])
                adi = ad(Tn0[:,:,1,arm]*expξ(xi[:,1,arm],q[i,1,arm]))
                for j = 2:n
                    A[(i-1)*6+1:6*i,(j-1)*6+1:6*j,arm] = adi*ad(Tn0[:,:,j,arm])
                    adi = adi*ad(Tn0[:,:,j,arm]*expξ(xi[:,j,arm],q[i,j,arm]))
                end
                A[(i-1)*6+1:6*i,n*6+1:6*(n+1),arm] = adi*ad(Tn0[:,:,n+1,arm])

            end
        end
        x = zeros(3)
        for arm = 1:2
            x = (A[:,:,arm]'A[:,:,arm] + λ*I)\A[:,:,arm]'y
            # Update the candidate nominal parameters
            s = arm == 1 ? 1 : -1 # Determine the sign of the correction
            for j = 1:n+1
                ξ = 0.5*s*x[(j-1)*6+1:6j] # Only apply half the update to each arm, verified to be essential
                Tn0cand[:,:,j,arm] = Tn0[:,:,j,arm]*expξ2(ξ,1)
            end

        end
        et[iter+1],er[iter+1] = evalErrorDual(xi,Tn0cand,q)
        println("Error: ",round(et[iter+1], digits=5), " λ: ", λ, " Norm dx: ", round(norm(x), digits=5))

        if norm(x) < 1e-14
            return Tn0,xi, et, er
        end
        if et[iter+1] > et[iter]
            et[iter+1] = et[iter]
            er[iter+1] = er[iter]
            λ *= 10
        else
            λ/=10
            xi = deepcopy(xini)
            Tn0 = deepcopy(Tn0cand)
        end

    end

    return Tn0, xi, et, er

end

"""
    xin, et, er = calibPOE(Xin,Ta,q; maxiter=50, λ = 10000.0)

Performs knematic calibration using the POE formulation of kinematics.
- `xin` Twists
- `et` Translational errors as function of iteration of algorithm
- `er` Rotational errors as function of iteration of algorithm
- `λ` Regularization parameter

See `runtests.jl` for simulated usage. Some functions to generate simulation data for calibration purposes are provided in $(@__FILE__)
"""
function calibPOE(Xin,Ta,q; maxiter=50, λ = 10000.0)

    xin     = copy(Xin)
    n       = size(xin,2)-2
    N       = size(q,1)
    y       = zeros(6N)
    A       = zeros(6N,6*(n+1))
    xini    = similar(xin)
    et      = zeros(maxiter+1)
    er      = zeros(maxiter+1)
    et[1],er[1]   = evalErrorPOE(xin,Ta,q)
    println("Error: ",et[1], " λ: ", λ)
    @assert size(Ta,3) == N

    for iter = 1:maxiter # do a few iterations of the calibration

        for i = 1:N
            Tfull = fkinePOE(xin,q[i,:]')
            # populate the matrices of the linear estimation problem
            y[(i-1)*6+1:6i] = twistcoords(logT(Ta[:,:,i]*trinv(Tfull)))
            A[(i-1)*6+1:6*i,:] = Ai(q[i,:]',xin)

        end

        x = (A'A + λ*I)\A'y
        # Update the candidate nominal parameters
        for j = 1:n+1
            ξ = x[(j-1)*6+1:6j]
            xini[:,j] = xin[:,j] + ξ
            #xini[:,j] = conformize(xini[:,j])
        end

        et[iter+1],er[iter+1] = evalErrorPOE(xini, Ta,q)
        println("Error: ",round(et[iter+1], digits=5), " λ: ", λ, " Norm dx: ", round(norm(x), digits=5))

        if norm(x) < 1e-10
            return xin, et, er
        end
        if et[iter+1] > et[iter]
            et[iter+1] = et[iter]
            λ *= 10
        else
            λ/=10
            xin = deepcopy(xini)
        end

    end

    return xin, et, er

end

function calibPOE_offsets_from_points(Xin,Q;maxiter=50, λ = 10000.0)
    xin     = copy(Xin)
    n       = size(xin,2)-2
    Ndatasets = size(Q,1)
    Npoints = sum(map(x-> size(x,1),Q))

    y       = zeros(6Npoints)
    A       = zeros(6Npoints,n-1)
    xini    = similar(xin)
    et      = zeros(maxiter+1)
    er      = zeros(maxiter+1)
    et[1]   = evalErrorPOE_offsets_from_points(xin,Q, zeros(n-1))
    println("Error: ",et[1], " λ: ", λ)
    δq      = zeros(n-1)

    for iter = 1:maxiter # do a few iterations of the calibration

        ii = 1

        for p = 1:Ndatasets
            q = Q[p]
            N = size(q,1)
            Tfull = cat(3,[fkinePOE(xin,q[i,:]'+[δq;0]) for i = 1:N]...)
            Tm = squeeze(mean(Tfull,3),3)
            toOrthoNormal!(Tm)
            for i = 1:N
                # populate the matrices of the linear estimation problem
                y[ii:ii+5] = twistcoords(logT(Tm*trinv(Tfull[:,:,i])))
                A[ii:ii+5,:] = xii(q[i,1:(n-1)]'+δq, xin[:,1:(n-1)])
                ii += 6
            end
        end

        x = (A'A + λ*I)\A'y
        # Update the candidate nominal parameters
        δqi = δq + x


        et[iter+1] = evalErrorPOE_offsets_from_points(xin, Q, δqi)
        println("$iter Error: ",round(et[iter+1], digits=5), " λ: ", λ, " Norm dx: ", round(norm(x), digits=5))

        if norm(x) < 1e-10
            return δq, et
        end
        if et[iter+1] > et[iter]
            et[iter+1] = et[iter]
            λ *= 10
        else
            # λ = max(λ/10,1e-8)
            δq = deepcopy(δqi)
        end

    end

    return δq, et

end

function evalError(xin,Tn0,Ta,q)
    N = size(Ta,3)
    et = 0.0
    er = 0.0
    for i = 1:N
        Tfull = fkineLPOE(Tn0,xin,q[i,:])
        et += norm(Ta[1:3,4,i]-Tfull[1:3,4])
        er += norm(skewcoords(logR(Ta[1:3,1:3,i]*Tfull[1:3,1:3]')))
    end
    return et/N, er/N
end

function evalErrorDual(xin,Tn0,q)
    N = size(Tn0,3)
    et = 0.0
    er = 0.0
    for i = 1:N
        T1 = fkineLPOE(Tn0[:,:,:,1],xin[:,:,1],q[i,:,1])
        T2 = fkineLPOE(Tn0[:,:,:,2],xin[:,:,2],q[i,:,2])
        et += norm(T2[1:3,4]-T1[1:3,4])
        er += norm(skewcoords(logR(T2[1:3,1:3]*T1[1:3,1:3]')))
    end
    return et/N, er/N
end

function evalErrorPOE(xin,Ta,q)
    xini = copy(xin)
    N = size(Ta,3)
    et = 0.0
    er = 0.0
    for i = 1:N
        Tfull = fkinePOE(xin,q[i,:])
        et += norm(Ta[1:3,4,i]-Tfull[1:3,4])
        er += norm(skewcoords(logR(Ta[1:3,1:3,i]*Tfull[1:3,1:3]')))
    end
    return et/N, er/N
end

function evalErrorPOE_offsets_from_points(xi,Q,dq)
    Ndatasets = size(Q,1)
    J = Array{Float64}(Ndatasets)
    for p = 1:Ndatasets
        q = Q[p]
        N = size(q,1)
        T = Array{typeof(dq[1])}(4,4,N)
        for i = 1:N
            T[:,:,i] = fkinePOE(xi,q[i,:] + [dq;0])
        end
        xyz = squeeze(T[1:3,4,:],2)
        mxyz = mean(xyz,2)
        Jpos = sum(abs((xyz.-mxyz)./mxyz))/N

        Rm = squeeze(mean(T[1:3,1:3,:],3),3)
        toOrthoNormal!(Rm)
        Jori = 0.0
        @inbounds for i = 1:N
            Jori += Rangle(T[1:3,1:3,i],Rm)
        end
        Jori /= Rangle(Rm)
        Jori /= N

        J[p] = Jpos + Jori
    end
    return sum(J)
end

"""
q, xin, T0, xinmod, Ta = simulateCalibration1(N)
Create simulated data
"""
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
        Ta[:,:,i] = fkinePOE(xin,q[i,:])
    end
    return q, xin, T0, xinmod, Ta
end

"""
q, xin, T0, xinmod, Ta = simulateCalibration_POE(N)
Create simulated data for use with `calibPOE`
"""
function simulateCalibration_POE(N)
    Random.seed!(1)
    n       = 6
    q       = 2π*rand(N,n)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, I4);
    xin     = DH2twistsPOE(dh)
    xinmod   = deepcopy(xin)
    for i = 1:n+1
        xinmod[:,i]  += [0.001randn(3); 0.01π/180*randn(3)]
        #conformize(xinmod[:,i])
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkinePOE(xin,q[i,:])
    end
    return q, xin, T0, xinmod, Ta
end

"""
q, xin, Tn0, Tn0mod, Ta = simulateCalibration_LPOE(N)
Create simulated data for use with `calibLPOE`
"""
function simulateCalibration_LPOE(N)
    Random.seed!(1)
    n       = 6
    q       = 2π*rand(N,n)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, I4);
    xin = DH2twistsLPOE(Tn0)
    Tn0mod = deepcopy(Tn0)
    for i = 1:n+1
        Tn0mod[1:3,4,i] += 0.1randn(3)
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkineLPOE(Tn0,xin,q[i,:])
        AAA,BBB,T = jacobian(q[i,:],dh, I4);
        @assert T ≈ Ta[:,:,i]
    end
    return q, xin, Tn0, Tn0mod, Ta
end

"""
q, xin, Tn0, Tn0mod, Ta = simulateCalibration_LPOE_dual(N)

Create simulated data for use with `calibLPOE_dual`
"""
function simulateCalibration_LPOE_dual(N)
    Random.seed!(1)
    n       = 6
    qt      = 2π*rand(N,n)
    q1      = qt + 0.01*π/180*randn(N,n)
    q2      = qt + 0.01*π/180*randn(N,n)
    q       = cat(q1,q2, dims=3)
    dh      = DH7600()
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, I4);
    xin = DH2twistsLPOE(Tn0)
    xin = cat(xin,xin, dims=3)
    Tn0mod = cat(deepcopy(Tn0),deepcopy(Tn0), dims=4)
    for i = 1:n+1
        Tn0mod[1:3,4,i,1] += 0.1randn(3)
        Tn0mod[1:3,4,i,2] += 0.1randn(3)
    end
    Ta  = zeros(4,4,N)
    for i = 1:N
        Ta[:,:,i] = fkineLPOE(Tn0,xin[:,:,1],qt[i,:])
    end

    return q, xin, Tn0, Tn0mod, Ta
end

#plot(et[et .!= 0])
