#Robotlib
This is a library of functions to help out in a robotics lab. At present stage, it contains functions for forward kinematics, jacobians, iterative inverse kinematics and for a few robotics related calibration problems. The library also contains a number of functions to convert from various orientation representations and other robotics related helper functions.

Install using

`Pkg.add("Robotlib")`

##Usage
```julia
fkine, ikine, jacobian = get_kinematic_functions("yumi") # Replace yumi for your robot model, as long as it's supported
data = orcalog2mat(pathopen, pathsave)
q = getData("robot_0.*posRawAbs", data, 1, removeNaN = false)
```

For YuMi, joint angles `q` must be converted to logical order using e.g. abb2logical!(q)
You must also consider the base transform of YuMi

The module imports the following files

```julia
include("POEutils.jl")
include("DH.jl")
include("kinematics.jl")
include("robotplot.jl")
include("Frames.jl")
include("read_log.jl")
include("csv2mat.jl")
```

##Exported functions
```julia
Rt2T, T2R, T2t, skewcoords, twistcoords, skew, skew4, expω, expξ, expξ2, expξ!, logT, logR
ad, adi, trinv, isrot, isse3, Rangle, conformize, DH2twistsPOE, DH2twistsLPOE, dh2Tn
toOrthoNormal!, toOrthoNormal, rpy2R, Quaternion, xyθ
fkinePOE, fkineLPOE, ikinePOE, dh2Tn, jacobianPOE, jacobianPOEikine, jacobian, get_kinematic_functions
plot_traj, plot_traj3, plot_traj_sub, plot3smart, plotsub
DH, DH2400, DHYuMi, DH7600, DHtest, abb2logical!, logical2abb!, abb2logical, logical2abb
csv2mat, orcalog2mat, getData, readmat
```

The module includes a submodule, Frames, which is aimed at replacing the Nikon K600 software. It supports creation of frames, simple projections, fitting of planes, lines etc. and has a number of plotting options. It must be separately imported with `using Robotlib.Frames`

The module includes a submodule, Calibration, which includes a number of calibration routines. It must be separately imported with `using Robotlib.Calibration`

## Kinematics
The library has functions for calculation of forward kinematics, inverse kinematics and jacobians. Several versions of all kinematics functions are provided; calculations can be made using either the DH-convention or the product of exponentials formulation. To support a new robot, create an object of the type `DH`, or provide a matrix with POE-style link twists, for use with the kinematic functions.
### Usage
```julia
dh = DH7600()
xi = DH2twistsPOE(dh)
T  = fkinePOE(xi, q)
```
or alternatively
```julia
dh = DH7600()
Jn, J0, T, Ti, trans = jacobian(q, dh)
```
many other options exits, check the kinematics.jl

#Frames
This module is aimed at assisting with the creation of frames for tracking using optical tracking systems. It supports projection of points and lines onto planes, creating frames from features and has some plotting functionality.

##Usage
```julia
path = "/home/fredrikb/work/flexifab/frames/"


add_frame_name!("SEAM","Weld seam frame")
add_frame_name!("TAB","Table frame")

T_RB_Tm = MAT.matread(path*"T_RB_T.mat")["T_RB_T"]
T_TF_TCPm = MAT.matread(path*"T_TF_TCP.mat")["T_TF_TCP"]
T_T_TABm = MAT.matread(path*"T_T_Table.mat")["T_T_Table"]

T_RB_T = Frame(T_RB_Tm,"RB","T")
T_S_D = Frame(T_TF_TCPm,"S","D")
T_T_TAB = Frame(T_T_TABm,"T","TAB")

cloud_seam = readcloud(path*"CloudSeam_edge.txt")
plane_seam = readplane(path*"PlaneSeam_edge.txt")
cloud_seam_projected = project(plane_seam,cloud_seam)
line_seam = fitline(cloud_seam_projected)

T_T_SEAM = framefromfeatures(("z+",line_seam),("y-",plane_seam),cloud_seam_projected[1],"SEAM")
T_RB_SEAM = T_RB_T*T_T_SEAM
T_RB_TAB = T_RB_T*T_T_TAB
T_TAB_SEAM = inv(T_T_TAB)*T_T_SEAM

MAT.matwrite(path*"T_TAB_SEAM.mat",["T_TAB_SEAM" => T_TAB_SEAM.T])
MAT.matwrite(path*"T_T_SEAM.mat",["T_T_SEAM" => T_T_SEAM.T])
MAT.matwrite(path*"T_RB_TAB.mat",["T_RB_TAB" => T_RB_TAB.T])
println("Wrote T_TAB_SEAM, T_T_SEAM, T_RB_TAB to files in $path")

cloud_seam_RB = T_RB_T*cloud_seam
cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
plane_seam_RB = T_RB_T*plane_seam
line_seam_RB = T_RB_T*line_seam



plotframe(Frame(eye(4),"RB","U"),200, label=true)

plotpoints(cloud_seam_RB,"bx")
plotpoints(cloud_seam_projected_RB,"rx")
plotline(line_seam_RB,"r",500,label="Line seam")
plotplane(plane_seam_RB,"b",200,label="Plane seam")
plotframe(T_RB_SEAM,200, label=true)
plotframe(T_RB_TAB,200, label=true)

PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")
PyPlot.axis("scaled")
```
