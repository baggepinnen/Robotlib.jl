using Plots
import .Frames # The dot is due to the fact that Frames is a submodule



@userplot TrajPlot
"""`trajplot(T,args...)` Plots a trajectory of T-matrices in a single plot"""
trajplot
@recipe function f(t::TrajPlot)
    title --> "Trajectory"
    layout --> 1
    x = t.args[1]
    squeeze(x[1:3,4,:],2)'
end


"""`plot_traj3(T, ls=".", plotFrame = 0)` Plots a trajectory of T-matrices in a 3D plot, with otional frames drawn of length `plotFrame`"""
function plot_traj3(T, ls=".", plotFrame = 0)
    plot(squeeze(T[1,4,:],(1,2)),squeeze(T[2,4,:],(1,2)),squeeze(T[3,4,:],(1,2)),ls)
    if plotFrame > 0
        for i = 1:size(T,3)
            Frames.plotframe(T[:,:,i],plotFrame)
        end
    end
end

@userplot Plot3Smart
"""`plot3smart(x,args...)` Makes a 3d plot of a matrix"""
plot3smart
@recipe function f(t::Plot3Smart)
    data = t.args[1]
    seriestype := :path3d
    (data[:,1][:],data[:,2][:],data[:,3][:])
end


# @recipe function f(::Type{Val{:mat3d}}, data)
#     @show typeof(data)
#     @show typeof(d[:x])
#     @show typeof(d[:y])
#     @show typeof(d[:z])
#     x := d[:y][:,1][:]
#     y := d[:y][:,2][:]
#     z := d[:y][:,3][:]
#     seriestype := :path3d
#     ()
# end
