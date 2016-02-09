import PyPlot
using Plots
import .Frames # The dot is due to the fact that Frames is a submodule


"""`plot_traj(T,args...)` Plots a trajectory of T-matrices in a single plot"""
plot_traj(T,args...;kw...) = plot(squeeze(T[1:3,4,:],2)',args...;kw...)


"""`plot_traj3(T, ls=".", plotFrame = 0)` Plots a trajectory of T-matrices in a 3D plot, with otional frames drawn of length `plotFrame`"""
function plot_traj3(T, ls=".", plotFrame = 0)
    PyPlot.plot3D(squeeze(T[1,4,:],(1,2)),squeeze(T[2,4,:],(1,2)),squeeze(T[3,4,:],(1,2)),ls)
    if plotFrame > 0
        for i = 1:size(T,3)
            Frames.plotframe(T[:,:,i],plotFrame)
        end
    end
end

"""`plot_traj_sub(T, args...)` plots a trajectory of T-matrices in three subplots"""
plot_traj_sub(T, args...;kw...) = subplot(squeeze(T[1:3,4,:],2)',args...; kw...)


"""`plot3smart(x,args...)` Makes a 3d plot of a matrix"""
plot3smart(x,args...) = PyPlot.plot3D(x[:,1],x[:,2],x[:,3],args...)


"""`plotsub(x, ls="-"; dir=:v, title="")` plots the columns of a matrix `x` in separate subplots"""
function plotsub(x, ls="-"; dir=:v, title="")
    (N,c) = size(x);

    if dir == :h || dir == "h"
        for i = 1:c
            PyPlot.subplot(1,c,i)
            PyPlot.plot(x[:,i],ls); PyPlot.hold(true); PyPlot.grid(true)

        end
        PyPlot.subplot(1,c,1)
        PyPlot.title(title)
    else
        for i = 1:c
            PyPlot.subplot(c,1,i)
            PyPlot.plot(x[:,i],ls); PyPlot.hold(true); PyPlot.grid(true)

        end
        PyPlot.subplot(c,1,1)
        PyPlot.title(title)
    end
end
