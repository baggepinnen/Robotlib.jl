[![Robotlib](http://pkg.julialang.org/badges/Robotlib_0.6.svg)](http://pkg.julialang.org/?pkg=Robotlib)
[![Build Status](https://travis-ci.org/baggepinnen/Robotlib.jl.svg?branch=master)](https://travis-ci.org/baggepinnen/Robotlib.jl)
[![codecov](https://codecov.io/gh/baggepinnen/Robotlib.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/baggepinnen/Robotlib.jl)


Now supporting Julia v0.7/v1.0.
Dependency MAT.jl is currently not building on Julia v1.0 and hence all file loading functionality and some tests have been deactivated.

# Robotlib
This is a library of functions to help out in a robotics lab. At present stage, it contains functions for forward kinematics, jacobians, iterative inverse kinematics and for a few robotics related calibration problems. The library also contains a number of functions to convert from various orientation representations and other robotics related helper functions.

Install using

`Pkg.add("Robotlib")`


## Usage
```julia
fkine, ikine, jacobian = get_kinematic_functions("yumi") # Replace yumi for your robot model, as long as it's supported
data = orcalog2mat(pathopen, pathsave) # Read data from a csv-file and save as binary file
q = getData("robot_0.*posRawAbs", data, 1, removeNaN = false) # Extract columns from data object using regex like syntax
```

For ABB YuMi, joint angles `q` must be converted to logical order using e.g. `abb2logical!(q)`
If you use the kinematic functions privided by `get_kinematic_functions`, the base transform is handled automatically. If you use the standard kinematic functions provided in Robotlib, you must also consider the base transform.

### Case study, calibrate force sensor
```julia
using Robotlib
using Robotlib.Calibration
using DSP # For filtfilt

# Define robot to use, in this case YuMi
dh = DHYuMi()
fkine, ikine, jacobian = get_kinematic_functions("yumi")

# Define paths to log file and where to store converted binary file for faster reading
pathopen = "/work/fredrikb/extRosetta/frida_gravity_2.txt"
pathsave = "/tmp/fredrikb/log.mat"

# Get data from the logfile
data    = orcalog2mat(pathopen, pathsave) # This line is only needed during the first run
data    = readmat(pathsave)
ds      = 1 # Downsampling factor
q       = getData("robot_1.*posRawAbs", data, ds) # Data vectors to retrieve are specified with regex style
q̇       = getData("robot_1.*velFlt", data, ds)
τ       = getData("robot_1.*trqRaw", data, ds)
f       = getData("force", data, ds)

# Convert joint data from ABB order to logical order
abb2logical!(q)
abb2logical!(q̇)
abb2logical!(τ)

# Apply gear ratio transformation
q = q*dh.GR'
q̇ = q̇*dh.GR'
τ = τ*inv(dh.GR')

# Filter velocities to get accelerations
q̈ = filtfilt(ones(50),[50.],centralDiff(q̇))

# plot(abs([q̇, q̈]))

# Sort out data with low acceleration
lowAcc  = all(abs.(q̈) .< 3e-4,2)
q       = q[lowAcc,:]
q̇       = q̇[lowAcc,:]
τ       = τ[lowAcc,:]
f       = f[lowAcc,:]
N       = size(q,1)


# Apply forward kinematics to get end-effector poses
T  = cat([fkine(q[i,:]) for i = 1:N]..., dims=3)

trajplot(T) # Plots a trajectory of R4x4 transformation matrices

# Perform the force sensor calibration and plot the errors
Rf,m,offset     = Robotlib.Calibration.calibForce(T,f,0.2205,offset=true)
err = hcat([Rf*f[i,1:3] + offset - T[1:3,1:3,i]'*[0, 0, m*-9.82] for i = 1:N]...)'
plot(f[:,1:3],lab="Force")
plot!(err,l=:dash,lab="Error")
println("Error: ", round(rms(err), digits=4))
```

## The module imports the following files

```julia
include("POEutils.jl")
include("DH.jl")
include("kinematics.jl")
include("robotplot.jl")
include("Frames.jl")
include("read_log.jl")
include("csv2mat.jl")
```

## Exported functions
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
The library has functions for calculation of forward kinematics, inverse kinematics and jacobians. Several versions of all kinematics functions are provided; calculations can be made using either the DH-convention or the (local) product of exponentials formulation. To support a new robot, create an object of the type `DH`, or provide a matrix with POE-style link twists, for use with the kinematic functions.
### Usage
```julia
dh = DH7600() # ABB Irb 7600
xi = DH2twistsPOE(dh)
T  = fkinePOE(xi, q)
```
or alternatively
```julia
dh = DH7600()
Jn, J0, T, Ti, trans = jacobian(q, dh)
```
many other options exits, check [kinematics.jl](src/kinematics.jl)

# Frames
This module is aimed at assisting with the creation of frames for tracking using optical tracking systems. It supports projection of points and lines onto planes, creating frames from features and has some plotting functionality.

## Usage
This is an example of how data can be loaded from files and how different geometrical objects can be fitted to data, projected onto other objects etc.
```julia
using Frames
import MAT
function setupframes(path)
	path = Pkg.dir("Robotlib","src","applications","frames/")

	# Add frame names to the dictionary
	add_frame_name!("SEAM","Weld seam frame")
	add_frame_name!("TAB","Table frame")

	# Read matrices from file
	T_RB_Tm = MAT.matread(path*"T_RB_T.mat")["T_RB_T"]
	T_TF_TCPm = MAT.matread(path*"T_TF_TCP.mat")["T_TF_TCP"]
	T_T_TABm = MAT.matread(path*"T_T_Table.mat")["T_T_Table"]

	# Create frames from matrices
	T_RB_T = Frame(T_RB_Tm,"RB","T")
	T_S_D = Frame(T_TF_TCPm,"S","D")
	T_T_TAB = Frame(T_T_TABm,"T","TAB")

	# Read point clouds generated by nikon software from file
	cloud_seam = readcloud(path*"CloudSeam_edge.txt")
	plane_seam = readplane(path*"PlaneSeam_edge.txt")

	# Project points onto plane and fit a line
	cloud_seam_projected = project(plane_seam,cloud_seam)
	line_seam = fitline(cloud_seam_projected)

	# Create a frame from the measured features
	T_T_SEAM = framefromfeatures(("z+",line_seam),("y-",plane_seam),cloud_seam_projected[1],"SEAM")
	T_RB_SEAM = T_RB_T*T_T_SEAM
	T_RB_TAB = T_RB_T*T_T_TAB
	T_TAB_SEAM = inv(T_T_TAB)*T_T_SEAM


	cloud_seam_RB = T_RB_T*cloud_seam
	cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
	plane_seam_RB = T_RB_T*plane_seam
	line_seam_RB = T_RB_T*line_seam

	# Plot results
	plot(Frame(I4,"RB","U"), 200)
	plot!(cloud_seam_RB, c=:blue)
	plot!(cloud_seam_projected_RB, c=:red)
	plot!(line_seam_RB, 500, label="Line seam")
	plot!(plane_seam_RB, 200, label="Plane seam")
	plot!(T_RB_SEAM, 200, label="T_RB_SEAM")
	plot!(T_RB_TAB, 200, label="T_RB_TAB")

	xlabel!("x")
	ylabel!("y")
	# zlabel!("z")

    # Write frames to file
    MAT.matwrite(path*"T_TAB_SEAM.mat",["T_TAB_SEAM" => T_TAB_SEAM.T])
    MAT.matwrite(path*"T_T_SEAM.mat",["T_T_SEAM" => T_T_SEAM.T])
    MAT.matwrite(path*"T_RB_TAB.mat",["T_RB_TAB" => T_RB_TAB.T])
    println("Wrote T_TAB_SEAM, T_T_SEAM, T_RB_TAB to files in $path")
end
```


# Citing
This software package was developed for the following thesis
[Bagge Carlson, F.](https://www.control.lth.se/staff/fredrik-bagge-carlson/), ["Machine Learning and System Identification for Estimation in Physical Systems"](https://lup.lub.lu.se/search/publication/ffb8dc85-ce12-4f75-8f2b-0881e492f6c0) (PhD Thesis 2018).
```bibtex
@thesis{bagge2018,
  title        = {Machine Learning and System Identification for Estimation in Physical Systems},
  author       = {Bagge Carlson, Fredrik},
  keyword      = {Machine Learning,System Identification,Robotics,Spectral estimation,Calibration,State estimation},
  month        = {12},
  type         = {PhD Thesis},
  number       = {TFRT-1122},
  institution  = {Dept. Automatic Control, Lund University, Sweden},
  year         = {2018},
  url          = {},
}
```
The algorithm `calibNAXP` was presented in
```bibtex
@inproceedings{bagge2015calibration,
  title        = {Six {DOF} eye-to-hand calibration from {2D} measurements using planar constraints},
  author       = {Bagge Carlson, Fredrik and Johansson, Rolf and Robertsson, Anders},
  booktitle    = {International Conference on Intelligent Robots and Systems (IROS)},
  year         = {2015},
  organization = {IEEE}
}
```
The friction model `frictionRBFN` was presented in
```bibtex
@inproceedings{bagge2015friction,
  title        = {Modeling and identification of position and temperature dependent friction phenomena without temperature sensing},
  author       = {Bagge Carlson, Fredrik and Robertsson, Anders and Johansson, Rolf},
  booktitle    = {International Conference on Intelligent Robots and Systems (IROS)},
  year         = {2015},
  organization = {IEEE}
}
```
