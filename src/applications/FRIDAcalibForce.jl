using Robotlib
using Robotlib.Calibration
using DSP # For filtfilt
using MAT

function demo_calibforce(filename = "")
    info("Running force calibration demo. This file uses the logfile gravitylog.mat and calculates the calibration parameters for a wrist mounted force sensor.")
    h          = 0.004032;
    dh         = DHYuMi()

    if if filename != "" # Use this option if you have a textbased logfile
        pathopen   = filename
        pathsave   = "log.mat"
        data       = orcalog2mat(pathopen, pathsave)
        # data       = readmat(pathsave)
        ds         = 1
        q          = getData("robot_1.*posRawAbs", data, ds)
        q̇          = getData("robot_1.*velFlt", data, ds)
        τ          = getData("robot_1.*trqRaw", data, ds)
        f          = getData("force", data, ds)
        # f[:,[2,5]] *= -1

        abb2logical!(q)
        abb2logical!(q̇)
        abb2logical!(τ)
        q = q*dh.GR'
        q̇ = q̇*dh.GR'
        τ = τ*inv(dh.GR')

        q̈ = filtfilt(ones(50),[50.],centralDiff(q̇))

        # plot(abs([q̇, q̈]))

        lowAcc = all(abs(q̈) .< 3e-4,2)[:]
        everyN = falses(size(q,1))
        everyN[1:10:end] = true
        q      = q[lowAcc & everyN,:]
        q̇      = q̇[lowAcc & everyN,:]
        τ      = τ[lowAcc & everyN,:]
        f      = f[lowAcc & everyN,:]
    else # Use this option to use the provided demofile
        data = readmat("data/gravitylog.mat")
        q = data["q"]
        q̇ = data["qd"]
        τ = data["tau"]
        f = data["f"]
    end

    N = size(q,1)
    fkine, ikine, jacobian = get_kinematic_functions("yumi")
    T = Array(Float64,4,4,N)
    for i = 1:N
        T[:,:,i]  = fkine(q[i,:])
    end

    # plot_traj(T)
    info("Calibrating force sensor")
    Rf,m,offset     = Robotlib.Calibration.calibForce(T,f,0.2205,offset=true)
    err = cat(2,[Rf*f[i,1:3] + offset - T[1:3,1:3,i]'*[0, 0, -9.82m] for i = 1:N]...)'
    plot(f[:,1:3],lab="Force")
    plot!(err,l=:dash,lab="Error")
    println("Error: ", round(rms(err),4))
    @show Rf
    @show m
    @show offset
    info("Done")

    # The code below can be used to save the result to a .mat-file
    # using MAT
    # matwrite("Rf.mat", Dict(
    #   	"Rf" => Rf,
    #     "mf" => m,
    #     "offset" => offset
    #   ))

    Rf,m,offset
end
