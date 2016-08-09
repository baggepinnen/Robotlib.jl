"""
This functions implements the algorithm from the paper
"Six DOF eye-to-hand calibration from 2D measurements using planar constraints"
which solves the problem `Prb = N'A(X*Ps)`
If result is bad, check if you send data in correct form\n
`POSES` is always the position of the tool frame from the robot FK\n
`points_S` are measurements from a line laser scanner. Given
in the sensor frame. This script assumes that the laser plane is the
sensor XY plane with y-axis pointing outwards from the sensor.\n
`lines_S` are vectors corresponding to the direction of the measured laser line
planes is a vector of indices corresponding to which plane a mesurement
comes from. It must be same length as `points_S`.\n
`N_calibs` determines the number of iterations\n
@inproceedings{carlson2015six,
  title={Six DOF eye-to-hand calibration from 2D measurements using planar constraints},
  author={Carlson, Fredrik Bagge and Johansson, Rolf and Robertsson, Anders},
  booktitle={Intelligent Robots and Systems (IROS), 2015 IEEE/RSJ International Conference on},
  pages={3628--3632},
  year={2015},
  organization={IEEE}
}
"""
function calibNAXP(points_S, lines_S, POSES, T_TF_S, planes,  N_calibs; doplot=false)
    N_planes = maximum(planes)
    RMScalibs = zeros(N_calibs)
    ALLcalibs = zeros(N_calibs,N_planes)
    for c = 1:N_calibs
        N_poses = size(POSES,3)
        # Convert points to RB cordinate system
        points = zeros(4,N_poses)
        lines = zeros(3,N_poses)
        for i = 1:N_poses
            T_RB_S = POSES[:,:,i] * T_TF_S
            points[:,i] = T_RB_S*[points_S[:,i]; 1]
            lines[:,i] = T2R(T_RB_S)*lines_S[:,i]
        end

        # Find plane centers, normals and distances
        mu_RB = zeros(3,N_poses)
        N_RB = mu_RB
        normals = zeros(N_planes,3)
        for j = 1:N_planes
            ind = planes .== j
            mu_t = mean(points[1:3,ind],2)
            mu_RB[:,ind] = repmat(mu_t,1,sum(ind))
            D,V = eig(cov(points[1:3,ind]'))
            N_RBt = V[:,1]
            if vecdot(mu_t,N_RBt) < 0
                N_RBt = -1*N_RBt
            end
            N_RBt = N_RBt*(N_RBt'mu_t)
            N_RB[:,ind] = repmat(N_RBt,1,sum(ind))
            normals[j,:] = N_RBt
        end

        # Estimation
        A = zeros(2*N_poses,9)
        y = zeros(2*N_poses,1)
        for i = 1:N_poses
            Nt = N_RB[:,i]
            Ra = POSES[1:3,1:3,i]
            Ta = POSES[1:3,4,i]
            Pt = points_S[:,i]
            y[i] = norm(N_RB[:,i])^2-vecdot(Nt,Ta)
            A[i,:] = (reshape(repmat(Nt'Ra,3,1)',9,1).*[reshape(repmat(Pt[1:2],1,3)',6,1);1;1;1])' #TODO: rewrite
            if true
                Pt = points_S[:,i] + 0.1*1.01^c*randn()*lines_S[:,i]
                y[i+N_poses] = norm(N_RB[:,i])^2-vecdot(Nt,Ta)
                A[i+N_poses,:] = (reshape(repmat(Nt'Ra,3,1)',9,1) .*[reshape(repmat(Pt[1:2],1,3)',6,1);1;1;1])'
            end

        end
        w = A\y
        er = y-A*w

        H = reshape(w,3,3)
        Rx = H[1:3,1]
        Ry = H[1:3,2]
        Rz = cross(Rx,Ry) # These lines must be changed if laser plane is changed
        Rnn = [Rx Ry Rz]
        @show Rnn
        @edb R = toOrthoNormal(Rnn)
        t = H[1:3,3]

        w2 = R[1:3,1:2]
        w2 = w2[:]
        y2 = y - A[:,1:6]*w2
        t2 = A[:,7:9]\y2
        if c < 0 #N_calibs
            T_TF_S = Rt2T(R,t)
        else
            T_TF_S = Rt2T(R,t2)
        end
        RMSi = zeros(N_planes)
        for j = 1:N_planes
            ind = planes .== j
            ind = find(ind)
            RMSi[j] = sqrt(pointDiff(T_TF_S,POSES[:,:,ind],points_S[1:3,ind])[1])
        end
        RMScalibs[c] = mean(RMSi)
        ALLcalibs[c,:] = RMSi'
        #  any(abs(RMSi) > 1e-1) && warn("Points does not seem to lie on a plane")
    end
    if doplot
        plotPlanes(normals)
        plotLines(points,lines)
    end
    norms = 0 # TODO: not implemented
    T_TF_S, RMScalibs,ALLcalibs, norms
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

function pointDiff(T_TF_S, T_RB_TF,points_S)
    loops = 1
    if size(T_TF_S,1) != 4
        loops = size(T_TF_S,2)
    end
    d = zeros(loops)
    for p = 1:loops
        if size(T_TF_S,1) != 4
            T_TF_S = Rt2T(rpy2R(T_TF_S[4:6,p]'),T_TF_S[1:3,p])
        end
        points = zeros(4,size(points_S,2))
        for i = 1:size(T_RB_TF,3)
            T_RB_S = T_RB_TF[:,:,i]*T_TF_S
            points[:,i] = T_RB_S*[points_S[:,i];1]
        end
        points = points[1:3,:]
        U,S,V = svd(points')
        d[p] = S[3]
        # d = max(SCOREp(:,3))^2;
        # d = mean(abs(SCOREp(:,3)))^2;
    end
    d
end
