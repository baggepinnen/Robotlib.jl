using Robotlib, TotalLeastSquares, HomotopyContinuation, DynamicPolynomials
using Test, Random
RTR(R) = R'R

# using SumOfSquares, PolyJuMP, Ipopt
# m = SOSModel(with_optimizer(Ipopt.Optimizer))

N         = 1000
Rf        = Robotlib.toOrthoNormal(randn(3,3))
gf        = randn(3)
POSES     = cat([Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
F    = hcat([Rf'*(POSES[1:3,1:3,i]'*gf) for i = 1:N]...)' |> copy
# function calibForceHomotopy(POSES,F,m0=1)
N = size(POSES,3)
@polyvar w[1:9]
@polyvar g[1:3]

# @variable(m, w[1:9])
# @variable(m, g[1:3])


R = reshape(w,3,3)
I₃ = Robotlib.I3
A  = Array{eltype(F)}(undef, 3N, 9)
B  = Array{Any}(undef, 3N)
At = (F,i) -> [F[i,1]*I₃   F[i,2]*I₃    F[i,3]*I₃]

for i = 1:N
    RA               = POSES[1:3,1:3,i]'
    b                = -RA*g
    A[3(i-1)+1:3i,:] = At(F,i)
    B[3(i-1)+1:3i]   = b
end
@info("Condition number of Gram matrix: ", cond(A'A))
@polyvar λ[1:6]
J = sum((A*w-B).^2) + sum(1:3) do i
    λ[i]*(R[:,i]'R[:,i] - 1)
end + λ[4]*R[:,1]'R[:,2] +
λ[5]*R[:,1]'R[:,3] +
λ[6]*R[:,2]'R[:,3]

dJ = DynamicPolynomials.differentiate(J, [w;g;λ])
result = HomotopyContinuation.solve(dJ)

reals = realsolutions(result)
costfun = w->sum((A*vec(w[1:9])-[b(w[10:12]) for b in B]).^2)
# R = tomat(reals[minindex])

result, J = calibForceHomotopy(POSES,forces)
reals = realsolutions(result)
tomat(x) = reshape(x[1:9],3,3) |> copy
Rs = tomat.(reals)
mininices = sortperm(J.(reals))
R̂ = Rs[mininices[1]]
ĝ = reals[mininices[1]][10:12]

# result, costfun, R, gf
# end



# set_start_value.(R, I₃)
# S = @set R'R .== I₃
# @constraint(m, R'R .== I₃, domain=S)
# @variable(m, γ)
# f = sum((A*w-B).^2)
# @constraint(m, f <= γ)
# # @objective(m, Min, γ)
# @objective(m, Min, f)
#
# optimize!(m)
# println(objective_value(m))
