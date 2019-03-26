using Robotlib, TotalLeastSquares, HomotopyContinuation, DynamicPolynomials
using Test, Random
RTR(R) = R'R

# using SumOfSquares, PolyJuMP, Ipopt
# m = SOSModel(with_optimizer(Ipopt.Optimizer))

# N         = 1000
# Rf        = Robotlib.toOrthoNormal(randn(3,3))
# mf        = 1
# gf        = [0,0,1] + 0.3randn(3)
# POSES     = cat([Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
# forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*gf) for i = 1:N]...)' |> copy
#
#
# function calibForceHomotopy(POSES,F,m0=1)
#     N = size(POSES,3)
#     @polyvar w[1:9]
#     @polyvar g[1:3]
#     @polyvar λ[1:6]
#     I = Robotlib.I3
#     A  = Array{eltype(F)}(undef, 3N, 9)
#     B  = Array{Any}(undef, 3N)
#     At = (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I]
#     for i = 1:N
#         RA               = POSES[1:3,1:3,i]'
#         A[3(i-1)+1:3i,:] = At(F,i)
#         B[3(i-1)+1:3i]   = -RA*g
#     end
#     @info("Condition number of Gram matrix: ", cond(A'A))
#     R = reshape(w,3,3)
#     J = sum((A*w-B).^2) + sum(1:3) do i
#         λ[i]*(R[:,i]'R[:,i] - 1)
#     end + λ[4]*R[:,1]'R[:,2] +
#     λ[5]*R[:,1]'R[:,3] +
#     λ[6]*R[:,2]'R[:,3]
#
#     dJ     = DynamicPolynomials.differentiate(J, [w;g;λ])
#     result = solve(dJ)
#     reals  = realsolutions(result)
#     costfun = w-> begin
#     sum((A*vec(w[1:9])-[b(w[10:12]) for b in B]).^2)
# end
# # R = tomat(reals[minindex])
# result, costfun
# end
#
#
#
# result, J = calibForceHomotopy(POSES,forces)
# reals = realsolutions(result)
# tomat(x) = reshape(x[1:9],3,3) |> copy
#
# Rs = tomat.(reals)
# RTR(R) = R'R
# minval, minindex = findmin(J.(reals))
# R = Rs[minindex]




function calibForceIterative(POSES,F,g; trace=false)
    N = size(POSES,3)
    I = Robotlib.I3
    A  = Array{eltype(F)}(undef, 3N, 3)
    B  = Array{eltype(F)}(undef, 3N)
    local m
    trace && (Rg = [])
    for iter = 1:6
        Rf, m = Robotlib.Calibration.calibForce(POSES,F,g; offset=false, verbose=false)
        for i = 1:N
            A[3(i-1)+1:3i,:] = Rf'POSES[:,:,i]'
            B[3(i-1)+1:3i]   = F[i,:]
        end
        g = A\B
        trace && push!(Rg, (Rf, g))
    end
    trace && (return Rf,g,m,Rg)
    Rf,g,m
end


##
σ = 10
traces = map(1:200) do mc
    N         = 1000
    Rf        = Robotlib.toOrthoNormal(randn(3,3))
    gf        = 100randn(3)
    POSES     = cat([Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
    forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*gf) for i = 1:N]...)' |> copy
    forces .+= σ*randn()

    ##

    @time R,g,m, Rg = calibForceIterative(POSES,forces,randn(3), trace=true)
    Rf,gf,Rg
end

errors = map(traces) do (Rf, gf, Rg)
    R,g = Rg[end]
    Rerr = Rangle(R,Rf, true)
    gerr = norm(g-gf) #< 1e-10 + √3*3σ/sqrt(N)
    Rerr,gerr
end

histogram([getindex.(errors, 1) getindex.(errors, 2)], layout=2)
##
errortraces = map(traces) do (Rf, gf, Rg)

    Rerr = [Rangle(R,Rf, true) for (R,g) in Rg]
    gerr = [norm(g-gf)/norm(gf) for (R,g) in Rg] #< 1e-10 + √3*3σ/sqrt(N)
    gaerr = [acos(g⋅gf/norm(gf)/norm(g))*180/pi for (R,g) in Rg]
    Rerr,gerr,gaerr
end

default(size=(800,600), grid=true, linealpha=0.3, linecolor=:black)
scales = (yscale=:log10, xscale=:log10, legend=false)
plot(getindex.(errortraces, 1); layout=(3,1), subplot=1, scales..., title="\$R\$ angle [deg]")
plot!(getindex.(errortraces, 2); subplot=2, scales..., title="\$g\$ relative error")
plot!(getindex.(errortraces, 3); subplot=3, scales..., title="\$g\$ angle [deg]")
# hline!([√3*3σ/sqrt(N)], l=(:black,:dash, 3), sp=2)
