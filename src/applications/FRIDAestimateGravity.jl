using Robotlib
# using Robotlib.Calibration
using DSP # For filtfilt
include("../dynamics.jl")
include("gravityFridaLS.jl")
include("../gravityFridaLS2.jl")
# import JuMP
# using Ipopt
h = 0.004032;


dh = DHYuMi()
# xi = DH2twistsPOE(dh)

pathopen = "/work/fredrikb/extRosetta/frida_gravity_2.txt"
pathsave = "/work/fredrikb/log.mat"

# pathopen = "frida_gravity_3.txt"
# pathsave = "log.mat"

# data    = orcalog2mat(pathopen, pathsave)
data    = readmat(pathsave);
ds      = 1
q       = getData("robot_1.*posRawAbs", data, ds);
q̇       = getData("irb2ext_robot_1.*velRef", data, ds);
τ       = getData("robot_1.*trqRaw", data, ds);
f       = getData("force", data, ds);


abb2logical!(q);
abb2logical!(q̇);
abb2logical!(τ);
q = q*dh.GR';
q̇ = q̇*dh.GR';
τ = τ*inv(dh.GR');

q̈ = filtfilt(ones(50),[50.],smartDiff(q̇));
q̇ = filtfilt(ones(10),[10.],q̇);
# plot(abs([q̇, q̈]))

lowAcc  = all(abs(q̈) .< 3e-4,2);
q       = q[lowAcc,:];
q̇       = q̇[lowAcc,:];
τ       = τ[lowAcc,:];
f       = f[lowAcc,:];
N       = size(q,1)

baseAnglesLeft  = [-0.63 , 0.95 , -0.18]
Rbase           = rpy2R(baseAnglesLeft,"xyz")
Tbase           = eye(4)
Tbase[1:3,1:3]  = Rbase
# T  = cat(3,[Tbase*fkinePOE(xi,q[i,:]') for i = 1:N]...);

# include("../dynamics.jl")

N = size(q,1)
np = 28+28
function getRegressor(q,q̇)
    N = size(q,1)
    A = Array(Float64,7N,np)
    ii = 1
    for i = 1:N
        A[ii:ii+6,:] = [gravityFridaLS(q[i,:])  diagm(q̇[i,:][:] .> 0)   diagm(q̇[i,:][:] .< 0) diagm(q̇[i,:][:] .> 0)*diagm(q̇[i,:][:])  diagm(q̇[i,:][:] .< 0)*diagm(q̇[i,:][:])]
        ii+= 7
    end
    return A
end

@time A = getRegressor(q,q̇);

still = vec(abs(q̇)') .< 2e-4
still[1:3:end] = true
still[2:3:end] = true

estA(A) = A[!still,:]
esty(y) = vec(y')[!still]


k = estA(A)\esty(τ)
# k = [estA(A); 0.1eye(np)]\[esty(τ);zeros(np)]
τhat = reshape(A*k,7,N)'
err = τ-τhat
println("Error: ", rms(esty(err)))



# println("Estimating a position dependent friction model for different regularization parameters")
# include("../posDepFric.jl")
# include("/local/home/fredrikb/SystemIdentification/src/least_squares.jl")
# n_basis = [10, 4]
# nppdf = prod(n_basis)
# centers = getCenters(n_basis, q, q̇)
# @time Apdf = frictionRBFN(q, q̇, centers, normalized=true);
# trials = 30
# λ = logspace(-15,6,trials)
# @time kpdf,errpdf,normk = dsvd(estA(Apdf),esty(err),λ)
# for i = 1:trials
#     τpdf = reshape(Apdf*kpdf[:,i],7,N)';
#     errpdf[i] = rms(esty(τ-τhat-τpdf))
# end
# Lcurve(errpdf, normk, λ)
# τpdf = reshape(Apdf*kpdf[:,end-8],7,N)';

plot(τ,layout=7)
plot!(τhat,c=:red,layout=7)
# subplot!(τhat + τpdf,c=:green,nr=7)
