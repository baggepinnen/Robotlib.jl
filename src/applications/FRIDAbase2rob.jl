using Robotlib
h = 0.004032;

info("This file calculates the rotation between the base and world coordinate systems.")

dh = DHYuMi()
xi = DH2twistsPOE(dh)

pathopen = "/work/fredrikb/extRosetta/base2rob2xyz.csv"
pathsave = "log.mat"

data    = orcalog2mat(pathopen, pathsave)
data    = readmat(pathsave)
ds      = 1
q       = getData("robot_1.*posRawAbs", data, ds)

# f[:,[2,5]] *= -1

abb2logical!(q)
q  = q*dh.GR'


baseAnglesLeft = [-0.63 , 0.95 , -0.18]
Rbase1         = rpy2R(baseAnglesLeft,"xyz")

N      = size(q,1)
T      = cat(3,[fkinePOE(xi,q[i,:]') for i = 1:N]...);
v1     = T[1:3,4,1500] - T[1:3,4,500]
v2     = T[1:3,4,2200] - T[1:3,4,1500]
v3     = T[1:3,4,3000] - T[1:3,4,2200]

v1    /= norm(v1)
v2    /= norm(v2)
v3    /= norm(v3)
Rbase2 = [v1 v2 v3]'

Rangle(Rbase1,Rbase2,true)
