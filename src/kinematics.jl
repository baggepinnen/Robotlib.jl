"""`Jn, J0, T, Ti, Tn = jacobian(qin, DH::DH, tool=eye(4))`
Calculates the jacobian in base (J0) and tool frame (Jn) as well as the forward kinematics (T),
given the joint angles and the DH-parameters"""
function jacobian{P}(q::AbstractVecOrMat{P}, DH::DH, tool=eye(4))
    n_joints    = size(q,1)
    dhpar       = DH.dhpar
    # q           = qin + DH.offset # Moved to dh2Tn
    trans       = dh2Tn(DH,q,tool)

    # ## jacobn (tool)
    Jn = zeros(P,6,n_joints)
    U = tool #Alternatively transformation to TCP <-------------- obs
    for j=n_joints:-1:1
        U = trans[:,:,j] * U
        # revolute axis
        d = [-U[1,1]*U[2,4]+U[2,1]*U[1,4];
        -U[1,2]*U[2,4]+U[2,2]*U[1,4];
        -U[1,3]*U[2,4]+U[2,3]*U[1,4]]
        delta = U[3,1:3]	# nz oz az

        Jn[:,j] = [d; delta]
    end

    ## fkine
    Ti          = zeros(4,4,n_joints+1)
    Ti[:,:,1]   = trans[:,:,1]
    T           = trans[:,:,1]
    for i=2:n_joints
        T           = T * trans[:,:,i]
        Ti[:,:,i] = Ti[:,:, i-1] * trans[:,:,i]
    end
    Ti[:,:,end] = Ti[:,:,end-1]*tool
    T           = T*tool

    # jacob0 (base)
    R = T[1:3,1:3]
    J0 = [R zeros(3,3); zeros(3,3) R] * Jn

    return Jn, J0, T, Ti, trans
end

"""
`jacobianPOE(q, xi)` Returns The jacobian in 1) the base frame, 2) the tool frame. It can most likely be rewritten to be faster.
"""
function jacobianPOE{P}(q::AbstractVecOrMat{P}, xi)
    #Page 117 in Murray 94
    n_joints  = size(q,1)
    Jss       = zeros(P,6,n_joints)
    adint     = eye(6)
    T         = eye(4)
    for i = 1:n_joints
        Jss[:,i] = adint*xi[:,i]
        Ti      = expξ(xi[:,i],q[i])
        T       = T*Ti
        adint   = adint*ad(Ti)
    end
    T = T * expξ(xi[:,end-1],1)*expξ(xi[:,end],1) # T is now the full FK

    # Now we do a bit messing around to get the jacobian we are used to from extCtrl
    Jbb = adi(T)*Jss
    Jsb = [T[1:3,1:3] zeros(3,3); zeros(3,3) T[1:3,1:3]] * Jbb

    return Jsb, Jbb, T
end

function jacobianPOEb{P}(q::AbstractVecOrMat{P}, xi)
    #Page 117 in Murray 94
    n_joints  = size(q,1)
    Jbb       = zeros(P,6,n_joints)
    T         = expξ2(xi[:,end-1],1)*expξ2(xi[:,end],1)
    adint     = adi(T)
    Ti        = eye(4)
    for i = n_joints:-1:1
        expξ2!(Ti,xi[:,i],q[i])
        adint   = adint*adi(Ti)
        Jbb[:,i]= adint*xi[:,i]
        T       = Ti*T
    end

    # Now we do a bit messing around to get the jacobian we are used to from extCtrl
    #Jbb = adi(T)*Jss
    Jsb = [T[1:3,1:3] zeros(3,3); zeros(3,3) T[1:3,1:3]] * Jbb

    return Jsb, Jbb, T
end

"""
`jacobianPOE(q, xi, T)` Returns The jacobian in 1) the base frame, 2) the tool frame and returns the FK. It must be handed a precomputed end transform `Tend = expξ2(xi[:,end-1],1)*expξ2(xi[:,end],1)`
"""
@fastmath function jacobianPOEikine(q, xi,T)
    #Page 117 in Murray 94
    n_joints  = size(q,1)
    Jbb       = zeros(6,n_joints)
    adint     = adi(T)
    Ti        = eye(4)
    @inbounds for i = n_joints:-1:1
        expξ2!(Ti,xi[:,i],q[i])
        adint   = adint*adi(Ti)
        Jbb[:,i]= adint*xi[:,i]
        T       = Ti*T
    end

    # Now we do a bit messing around to get the jacobian we are used to from extCtrl
    Jsb = [T[1:3,1:3] zeros(3,3); zeros(3,3) T[1:3,1:3]] * Jbb

    return Jsb, Jbb, T
end



"""Computes a set of local transformation matrices given the DH-parameters of a robot. Can be sent an optional joint-angle vector and a tool"""
function dh2Tn!{P}(Tn, DH, q::VecOrMat{P}=zeros(size(DH.dhpar,1)), tool=eye(4))
    dhpar    = DH.dhpar
    n_joints = size(dhpar,1)
    q = q .+ DH.offset
    for j=1:n_joints
        Tn[:,:,j] = [    cos(q[j])   -sin(q[j])*cos(dhpar[j,1])    sin(q[j])*sin(dhpar[j,1])    dhpar[j,2]*cos(q[j]);
        sin(q[j])   cos(q[j])*cos(dhpar[j,1])     -cos(q[j])*sin(dhpar[j,1])   dhpar[j,2]*sin(q[j]);
        0                 sin(dhpar[j,1])                     cos(dhpar[j,1])            dhpar[j,4];
        0                 0                                   0                                  1 ]
    end
    Tn[:,:,end] = tool
    return Tn
end

function dh2Tn{P}(DH, q::VecOrMat{P}=zeros(size(DH.dhpar,1)), tool=eye(4))
    n_joints = size(DH.dhpar,1)
    Tn       = zeros(P,4,4,n_joints+1)
    return dh2Tn!(Tn, DH, q, tool)
end

"""`fkinePOE(xi0,q)` Forward kinematics using POE"""
function fkinePOE(xi0,q)
    Tfull = expξ(xi0[:,end-1],1)*expξ(xi0[:,end],1) # We first handle the tool transform, which does not include a joint variable
    for j = size(xi0,2)-2:-1:1
        ξ   = xi0[:,j]
        Tfull = expξ(ξ,q[j])*Tfull
    end
    return Tfull
end

"""`fkineLPOE(Tn0,xi,q)` Forward kinematics using LPOE"""
function fkineLPOE(Tn0,xi,q)
    T = eye(4)
    n = size(xi,2)-1
    for j = 1:n
        T = T*Tn0[:,:,j]*expξ(xi[:,j],q[j]) # turn all joints
    end
    T = T*Tn0[:,:,end]
    return T
end

"""
`function ikinePOE(xi,T,q0; maxiter=100, λ = 1e0, tol = 1e-12, verbose = false, adaptive = true)`
Iterative inverse kinematics
`λ`sets the start value of the regularizing parameter (kind of like inverse step size). Bigger value means slower but more robust convergence, if adaptive is set to true, λ is adapted to speed up convergence.
"""
function ikinePOE(xi,T,q0; maxiter=100, λ = 1e0, tol = 1e-12, verbose = false, adaptive = true)
    q = copy(q0)
    err = Array{eltype(q0)}(maxiter)
    Tend = expξ2(xi[:,end-1],1)*expξ2(xi[:,end],1)
    J,V,Tc = jacobianPOEikine(q,xi,Tend)
    err[1] = norm(twistcoords(logT(T*trinv(Tc))))
    newjac = true
    iter = 2
    while iter <= maxiter

        dq = (J'J + λ*I)\J'*twistcoords(logT(T*trinv(Tc)))
        qc = q + dq
        J,V,Tc = jacobianPOEikine(qc,xi,Tend)
        err[iter] = norm(twistcoords(logT(T*trinv(Tc))))
        if err[iter] > err[iter-1]
            if adaptive
                λ *= 8
            end
            err[iter] = err[iter-1]
            newjac = false
        else
            q = qc
            if adaptive
                λ /= 2
            end
            newjac = true
        end
        verbose && println("Error: ",round(err[iter],7), " λ: ", round(λ,10), " Norm dq: ", round(norm(dq),7))
        if err[iter] < tol || norm(dq) < tol
            err = err[1:iter]
            break
        end

        iter += 1
    end
    return q, err
end

"""
Run this function with a string representing the robot you want the kinematic functions for, e.g. `get_kinematic_functions("YuMi")`
Currently supports YuMi and ABB IRB7600
returns `fkine(q), ikine(T,q0), jacobian(q)`
TODO: implement YuMi_left, Yumi_right Tbase*ikinePOE
"""
function get_kinematic_functions(robot)
    robot = lowercase(robot)
    if contains(robot,"yumi") || contains(robot,"frida")
        dh = DHYuMi()
    elseif robot == "7600" || robot == "irb7600"
        dh = DH7600()
    end
    Tn0 = dh2Tn(dh)
    xi = DH2twistsPOE(dh)
    xiL = DH2twistsLPOE(dh)
    Z = zeros(3,3)
    if robot == "yumileft"
        baseAnglesLeft = Quaternion(0.82888, -0.31402, 0.40801, -0.2188)
        TbaseLeft = [rotationmatrix(baseAnglesLeft) [0.0476, 0.07, 0.4115]; 0 0 0 1]
        fkinef = q -> TbaseLeft*fkineLPOE(Tn0,xiL,q)
        jacobianf = q -> [TbaseLeft[1:3,1:3] Z;Z TbaseLeft[1:3,1:3]]*jacobianPOE(q,xi)[1]
        ikinef = (T,q0, maxiter=100, λ = 1e0, tol = 1e-12, verbose = false) -> ikinePOE(xi,trinv(Tbase)*T,q0,maxiter=maxiter, λ = λ, tol = tol, verbose = verbose)
    elseif robot == "yumiright"
        baseAnglesRight = Quaternion(0.82888, 0.31402, 0.40801, 0.2188)
        TbaseRight = [rotationmatrix(baseAnglesRight) [0.0476, -0.07, 0.4115]; 0 0 0 1]
        fkinef = q -> TbaseRight*fkineLPOE(Tn0,xiL,q)
        jacobianf = q -> [TbaseRight[1:3,1:3] Z;Z TbaseRight[1:3,1:3]]*jacobianPOE(q,xi)[1]
        ikinef = (T,q0, maxiter=100, λ = 1e0, tol = 1e-12, verbose = false) -> ikinePOE(xi,trinv(Tbase)*T,q0,maxiter=maxiter, λ = λ, tol = tol, verbose = verbose)
    else
        fkinef = q -> fkineLPOE(Tn0,xiL,q)
        jacobianf = q -> jacobianPOE(q,xi)[1]
        ikinef = (T,q0, maxiter=100, λ = 1e0, tol = 1e-12, verbose = false) -> ikinePOE(xi,T,q0,maxiter=maxiter, λ = λ, tol = tol, verbose = verbose)
    end
    return fkinef, ikinef, jacobianf
end


function testikine()
    dh = DH7600();
    q0 = randn(6)
    xi = DH2twistsPOE(dh)
    T0 = fkinePOE(xi,q0)
    Tt = deepcopy(T0)
    Tt[1:3,4] += 0.1*[1,1,1]
    Tt[1:3,1:3] *= expω(1*π/180*randn(3))
    Profile.clear()
    gc()
    @time q,err = ikinePOE(xi,Tt,q0,verbose=true)

    # pp = semilogy(err)
    # display(pp)
end


function testJacobian()
    include("../POEutils.jl")
    dh = DH7600();
    q = [1,1,0,0,0,0];
    q = randn(6)
    AAA,BBB,T0,Ti0,Tn0 = jacobian(zeros(6),dh, eye(4));
    Jn,J0,Tjac,Tijac,Tnjac = jacobian(q,dh, eye(4));

    xiL = DH2twistsLPOE(Tn0)
    xi = DH2twistsPOE(Tn0)
    Jsb, Jbb = jacobianPOE(q,xi)

    display(round(J0,3))
    display(round(Jsb,3))
end
