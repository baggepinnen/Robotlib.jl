using Robotlib, TotalLeastSquares, FillArrays
# using HomotopyContinuation, DynamicPolynomials
using Test, Random
RTR(R) = R'R
vecangle(g,gf) = acos(min(1,g⋅gf/norm(gf)/norm(g)))*180/pi


rcond(K) = /(svdvals(K)[[1,3]]...)
toR(r) = toOrthoNormal(reshape(real(r),3,3))
tor(r) = toR(r)[:]

function rqi(K,r,n=1)
    for i = 1:n
        r = (K-I)\r
        r = r/norm(r)
    end
    r
end

function eigenR(K)
    v = real(eigen(K).vectors[:,1])
    R = reshape(v,3,3)
    det(R) < 0 && (R .*= -1)
    toOrthoNormal(R)
end
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
    B = F'[:]
    local m
    trace && (Rg = [])
    for iter = 1:6
        Rf, m = Robotlib.Calibration.calibForce(POSES,F,g; offset=false, verbose=true)
        for i = 1:N
            A[3(i-1)+1:3i,:] = Rf'POSES[:,:,i]'
        end
        g = A\B
        trace && push!(Rg, (Rf, g))
    end
    trace && (return Rf,g,m,Rg)
    Rf,g,m
end

function calibForceIterative2(POSES,F,g; trace=false)
    N = size(POSES,3)
    I = Robotlib.I3
    A  = Array{eltype(F)}(undef, 3N, 3)
    for i = 1:N
        A[3(i-1)+1:3i,:] = POSES[:,:,i]'
    end
    A = factorize(A)
    local m
    trace && (Rg = [])
    for iter = 1:6
        Rf, m = Robotlib.Calibration.calibForce(POSES,F,g; offset=false, verbose=false)
        B = (Rf*F')[:]
        g = A\B
        trace && push!(Rg, (Rf, g))
    end
    trace && (return Rf,g,m,Rg)
    Rf,g,m
end


function calibForceIterative3(POSES,F; trace=false)
    N = size(POSES,3)
    I = Robotlib.I3
    Rf = I
    A  = Array{eltype(F)}(undef, 3N, 3)
    B  = Array{eltype(F)}(undef, 3N)
    local m, g
    trace && (Rg = [])
    for iter = 1:6
        for i = 1:N
            A[3(i-1)+1:3i,:] = Rf'POSES[:,:,i]'
            B[3(i-1)+1:3i]   = F[i,:]
        end
        g = A\B
        Rf, m = Robotlib.Calibration.calibForce(POSES,F,g; offset=false, verbose=false)
        trace && push!(Rg, (Rf, g))
    end
    trace && (return Rf,g,m,Rg)
    Rf,g,m
end
using FillArrays

function calibForceIterative4(POSES,forces; trace=false)
    N = size(POSES,3)
    I3 = Matrix{Float64}(I, 3, 3)
    r = vec(I3)
    # D = permutedims(POSES, (2,1,3))
    # D = reshape(D, 3N,3)
    D = reshape(POSES, 3,3N)'
    F = kron(forces, I3)
    FF = pinv(F)*D
    DD = pinv(D)*F
    K = FF*DD
    e = eigen(K)
    sort(e.values, by=abs,rev=true)
    trace && (Rg = [])
    for iter = 1:5
        r = K*r
        r = vec(toOrthoNormal(reshape(r,3,3)))
        # trace && push!(Rg, (Rf, g))
    end
    # trace && (return Rf,g,m,Rg)
    # Rf,g,m
    toOrthoNormal(reshape(r,3,3)), e
end

function calibForceEigen(POSES,forces)
    N = size(POSES,3)
    D = reshape(POSES, 3,3N)'
    F = kron(forces, Eye(3))
    K = F\D*(D\F)
    R̂ = eigenR(K)
    ĝ = D\F*vec(R̂)
    R̂,ĝ
end



##
σ = 1
traces = map(1:30) do mc
    N         = 100
    Rf        = Robotlib.toOrthoNormal(randn(3,3))
    gf        = 100randn(3)
    POSES     = cat([Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
    forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*gf) for i = 1:N]...)' |> copy
    forces .+= σ*randn(size(forces))
    R,g,m, Rg = calibForceIterative3(POSES,forces, trace=true)
    R,g = calibForceEigen(POSES,forces)
    Rf,gf,Rg, R,g
end
##
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
    gaerr = [vecangle(g,gf) for (R,g) in Rg]
    Rerr,gerr,gaerr
end

default(size=(800,600), grid=true, linealpha=0.3, linecolor=:black)
scales = (yscale=:identity, xscale=:identity, legend=false)
ylims = (ylims=(0,30),)
plot(getindex.(errortraces, 1); layout=(1,3), subplot=1, scales..., title="\$R\$ angle [deg]")
hline!([180], subplot=1, l=(:black,:dash))
plot!(getindex.(errortraces, 2); subplot=2, scales..., title="\$g\$ relative error")
plot!(getindex.(errortraces, 3); subplot=3, scales..., title="\$g\$ angle [deg]")

##
opts = (fillalpha=0.3,)
histogram([e[1][end] for e in errortraces]; layout=(3,1), subplot=1, title="Rangle", opts...)
histogram!([Rangle(t[1],t[4])*180/pi for t in traces]; subplot=1, opts...)

histogram!([e[2][end] for e in errortraces]; subplot=2, title="g relative", opts...)
histogram!([norm(t[2]-t[5])/norm(t[2]) for t in traces]; subplot=2, opts...)

histogram!([e[3][end] for e in errortraces]; subplot=3, title="g angle", opts...)
histogram!([vecangle(t[2],t[5]) for t in traces]; subplot=3, opts...)

# hline!([√3*3σ/sqrt(N)], l=(:black,:dash, 3), sp=2)

##
using Test
A = 10randn(100,3)
B = 10randn(100,3)

@test 1/opnorm(A) < opnorm(pinv(A))
# @show sort(abs.(eigvals(A*pinv(A)*B*pinv(B))), rev=true)

##
D = [toOrthoNormal(randn(3,3)) for _ in 1:1000]
B = [D' for D in D]
D = reduce(vcat,D)
B = reduce(vcat,B)
opnorm(D)
opnorm(B*pinv(D))
sort(eigvals(B*pinv(D)), by=abs, rev=true)
svdvals(B*pinv(D))




##
σ = 1
N         = 1000
Rf        = Robotlib.toOrthoNormal(randn(3,3))
gf        = randn(3)
gf        = gf/norm(gf)
A         = [Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]
POSES     = reduce((x,y)->cat(x,y, dims=3), A)
D         = reshape(POSES, 3,3N)'
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*gf) for i = 1:N]...)' |> copy
forces  .+= σ*randn()
F = kron(forces, Eye(3))

# FF = F*pinv(F)
# DD = D*pinv(D)

K = F\D*(D\F)
# K = pinv(F)*D*(pinv(D)*F)
# eigvals(K)
# svdvals(K)
# opnorm(K)
# rcond(K)
Rh = eigenR(K)
norm((D\F*vec(Rh)) - gf)
# Rangle(Rh,Rf)*180/pi

##
mapslices(eigen(K).vectors, dims=1) do v
    vecangle(real.(v),vec(Rf))
end

r = vec(Eye(3)) |> Vector
for i = 1:5
    global r
    # r = tor(K*r)
    r .= K*r
    r ./= norm(r)
    @show norm(toR(r) - Rf)
end
println()

R = reshape(real(eigen(K).vectors[:,1]),3,3)
R'Rf
