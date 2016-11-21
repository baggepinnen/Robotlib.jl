using Robotlib
using Robotlib.Calibration
using DSP # For filtfilt

h          = 0.004032;
dh         = DHYuMi()
xi         = DH2twistsPOE(dh)
# pathopen   = "data/frida_gravity.txt"
pathsave   = "log.mat"
# data       = orcalog2mat(pathopen, pathsave)
data       = readmat(pathsave)
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

lowAcc  = all(abs(q̈) .< 3e-4,2)[:]
q       = q[lowAcc,:]
q̇       = q̇[lowAcc,:]
τ       = τ[lowAcc,:]
f       = f[lowAcc,:]
N   = size(q,1)

baseAnglesLeft  = [-0.63 , 0.95 , -0.18]
Rbase           = rpy2R(baseAnglesLeft,"xyz")
Tbase           = eye(4)
Tbase[1:3,1:3]  = Rbase

fkine, ikine, jacobian = get_kinematic_functions("yumi")
T = Array(Float64,4,4,N)
for i = 1:N
    T[:,:,i]  = Tbase*fkinePOE(xi,q[i,:]')
end


# plot_traj(T)


Rf,m,offset     = Robotlib.Calibration.calibForce(T,f,0.2205,offset=true)
err = cat(2,[Rf*f[i,1:3]' + offset - T[1:3,1:3,i]'*[0, 0, m*-9.82] for i = 1:N]...)'
plot(f[:,1:3],lab="Force")
plot!(err,l=:dash,lab="Error")
println("Error: ", round(rms(err),4))



# using MAT
# matwrite("Rf.mat", Dict(
#   	"Rf" => Rf,
#     "mf" => m,
#     "offset" => offset
#   ))
