import .Frames # The dot is due to the fact that Frames is a submodule
using RecipesBase
# default(markersize=1)

@userplot TrajPlot
"""`trajplot(T,args...)` Plots a trajectory of T-matrices in a single plot"""
trajplot
@recipe function trajplot(t::TrajPlot)
    title --> "Trajectory"
    layout --> 1
    T = t.args[1]
    label --> ["x" "y" "z"]
    T[1:3,4,:]'
end


@userplot TrajPlot3
"""`plot_traj3(T, ls=".", plotFrame = 0)` Plots a trajectory of T-matrices in a 3D plot, with otional frames drawn of length `plotFrame`"""
trajplot3
@recipe function trajplot3(t::TrajPlot3, plotFrame = 0)
    T = t.args[1]
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    @series begin
        seriestype := :scatter3d
        (T[1,4,:],T[2,4,:],T[3,4,:])
    end
    if plotFrame > 0
        for i = 1:size(T,3)
            @series begin
                plot!(T[:,:,i],plotFrame)
            end
        end
    end
end
