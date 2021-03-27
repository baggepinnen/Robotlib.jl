"""
See documentation at https://github.com/baggepinnen/Robotlib.jl
"""
module Robotlib
using LinearAlgebra, Statistics, StaticArrays, SparseArrays
using TotalLeastSquares
using Quaternions
import Quaternions: Quaternion, rotationmatrix
export rotationmatrix
const I4 = SMatrix{4,4,Float64,16}(Matrix{Float64}(I, 4, 4))
const I3 = SMatrix{3,3,Float64,9}(Matrix{Float64}(I, 3, 3))
include("DH.jl")
include("utils.jl")
include("kinematics.jl")
include("Frames.jl") # Must come before robotplot
include("robotplot.jl")
include("read_log.jl")
include("posDepFric.jl")



include("calibLPOE.jl")
include("calib_force.jl")
include("calibNAXP.jl")
include("calibAXYB.jl")

export calibLPOE,
    calibLPOEdual, calib_force, calib_force_iterative, calib_force_eigen, calibNAXP, calibAXYB, leds_to_tf

using Requires
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        function plotLines(points, lines)
            P = size(points, 2)
            for i = 1:P
                p = points[1:3, i]
                l = lines[1:3, i]
                Plots.plot3d!(
                    eachcol([p' + 0.05 * l'; p' - 0.05 * l'])...,
                    c = :red,
                    lab = "",
                )
            end
            display(Plots.current())
        end

        function plotPlanes(normals)
            P = size(normals, 1)
            p = zeros(P, 4)
            Plots.plot()
            for i = 1:P
                n = normals[i, :]
                xdir = n + [1, 0, 0]
                ydir = normalize(cross(n, xdir))
                xdir = normalize(cross(ydir, n))
                p[:, 1] = n + xdir + ydir
                p[:, 2] = n + xdir - ydir
                p[:, 3] = n - xdir - ydir
                p[:, 4] = n - xdir + ydir

                Plots.surface!(p[1, :], p[2, :], p[3, :])
            end
            Plots.plot!(alpha = 0.5, zlims = (-0.3, 1))
        end
    end
end





export skewcoords,
    twistcoords,
    skew,
    skew4,
    expω,
    expξ,
    expξ2,
    expξ!,
    logT,
    logR,
    trinv,
    isrot,
    isse3,
    Rangle,
    conformize,
    DH2twistsPOE,
    DH2twistsLPOE,
    dh2Tn,
    orthonormal,
    orthonormal!,
    rpy2R,
    Quaternion,
    xyθ,
    smartdiff,
    R2rpy,
    traj2quat,
    centraldiff
export fkinePOE,
    fkineLPOE,
    ikinePOE,
    dh2Tn,
    jacobianPOE,
    jacobianPOEikine,
    jacobian,
    get_kinematic_functions
export trajplot, trajplot3, plot3smart
export DH,
    DH2400,
    DHYuMi,
    DH7600,
    DHtest,
    abb2logical!,
    logical2abb!,
    abb2logical,
    logical2abb
export csv2dict, csv2mat, orcalog2mat, getdata
export frictionRBFN, getcenters

dir(x...) = joinpath(dirname(pathof(Robotlib)), x...)


end