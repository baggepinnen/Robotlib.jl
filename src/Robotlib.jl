"""
This is a library of functions to help out in the Robotlab at LTH
The module imports the following files\n
include("utils.jl")\n
include("DH.jl")\n
include("kinematics.jl")\n
include("robotplot.jl")\n
include("Frames.jl")\n
include("read_log.jl")\n
include("csv2mat.jl")\n
include("posDepFric.jl")\n

The module includes a submodule, Frames, which is aimed at replacing the Nikon K600 software. It supports creation of frames, simple projections, fitting of planes, lines etc. and has a number of plotting options. It must be separately imported with `using Robotlib.Frames`

The module includes a submodule, Calibration, which includes a number of calibration routines. It must be separately imported with `using Robotlib.Calibration`

Usage:
fkine, ikine, jacobian = get_kinematic_functions("yumi")
data = orcalog2mat(pathopen, pathsave)
q = getData(\"robot_0.*posRawAbs\", data, 1, removeNaN = false)

For YuMi, joint angles `q` must be converted to logical order using e.g. abb2logical!(q)
You must also consider the base transform of YuMi
"""
module Robotlib
using LinearAlgebra, Statistics, StaticArrays, SparseArrays
const I4 = SMatrix{4, 4, Float64, 16}(Matrix{Float64}(I, 4, 4))
const I3 = SMatrix{3, 3, Float64, 9}(Matrix{Float64}(I, 3, 3))
include("DH.jl")
include("utils.jl")
include("kinematics.jl")
include("Frames.jl") # Must come before robotplot
include("robotplot.jl")
# https://github.com/JuliaIO/MAT.jl/issues/90
include("read_log.jl")# Re-enable when MAT.jl builds under Julia v1.0. Then add MAT to require
include("csv2mat.jl") # Re-enable when MAT.jl builds under Julia v1.0. Then add MAT to require
include("posDepFric.jl")

module Calibration
using TotalLeastSquares
using ..Robotlib
using ..Robotlib: T2R,Rt2T, T2R, T2t, skewcoords, twistcoords, I3, I4
using LinearAlgebra, Statistics, SparseArrays
import Robotlib.xii
import Robotlib.Ai
include("calibLPOE.jl")
include("calibForce.jl")
include("calibNAXP.jl")
export calibLPOE, calibLPOEdual, calibForce, calibPOE_offsets_from_points, calibNAXP

using Requires
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        function plotLines(points,lines)
            P = size(points,2)
            for i = 1:P
                p = points[1:3,i]
                l = lines[1:3,i]
                plot3smart!([p'+0.05*l'; p'-0.05*l'],c=:red, lab="")
            end
            display(Plots.current())
        end

        function plotPlanes(normals)
            P = size(normals,1)
            p = zeros(P,4)
            Plots.plot()
            for i = 1:P
                n      = normals[i,:]
                xdir   = n+[1, 0, 0]
                ydir   = normalize(cross(n,xdir))
                xdir   = normalize(cross(ydir,n))
                p[:,1] = n + xdir + ydir
                p[:,2] = n + xdir - ydir
                p[:,3] = n - xdir - ydir
                p[:,4] = n - xdir + ydir

                Plots.surface!(p[1,:],p[2,:],p[3,:])
            end
            Plots.plot!(alpha=0.5, zlims=(-0.3, 1))
        end
    end
end


end


export skewcoords, twistcoords, skew, skew4, expω, expξ, expξ2, expξ!, logT, logR, trinv, isrot, isse3, Rangle, conformize, DH2twistsPOE, DH2twistsLPOE, dh2Tn, toOrthoNormal, toOrthoNormal!, rpy2R, Quaternion, xyθ, smartDiff, R2rpy, traj2quat, centraldiff
export fkinePOE, fkineLPOE, ikinePOE, dh2Tn, jacobianPOE, jacobianPOEikine, jacobian, get_kinematic_functions
export trajplot, trajplot3, plot3smart
export DH, DH2400, DHYuMi, DH7600, DHtest, abb2logical!, logical2abb!, abb2logical, logical2abb
export csv2mat, orcalog2mat, getData, readmat
export frictionRBFN, getCenters

dir(x...) = joinpath(dirname(pathof(Robotlib)),x...)


precompile(get_kinematic_functions, (String,))
precompile(fkineLPOE, (Array{Float64,3},Array{Float64,2},Array{Float64,1}))
precompile(jacobianPOE, (Array{Float64,1},Array{Float64,2}))



end
