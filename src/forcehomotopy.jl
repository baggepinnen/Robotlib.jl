using Robotlib, TotalLeastSquares, HomotopyContinuation, DynamicPolynomials
using Test, Random

N         = 1000
Rf        = Robotlib.toOrthoNormal(randn(3,3))
mf        = 1
g = randn(3)
POSES     = cat([Robotlib.toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*g) for i = 1:N]...)' |> copy


function calibForceHomotopy(POSES,F,m0=1)

    N = size(POSES,3)

    @polyvar w[1:9]
    @polyvar g[1:3]
    @polyvar λ[1:6]
    I = Robotlib.I3
    A  = Array{eltype(F)}(undef, 3N, 9)
    B  = Array{Any}(undef, 3N)
    At = (F,i) -> [F[i,1]*I   F[i,2]*I    F[i,3]*I]

    for i = 1:N
        RA               = POSES[1:3,1:3,i]'
        b                = -RA*g
        A[3(i-1)+1:3i,:] = At(F,i)
        B[3(i-1)+1:3i]   = b
    end
    @info("Condition number of Gram matrix: ", cond(A'A))
    R = reshape(w,3,3)
    J = sum((A*w-B).^2) + sum(1:3) do i
        λ[i]*(R[:,i]'R[:,i] - 1)
    end + λ[4]*R[:,1]'R[:,2] +
    λ[5]*R[:,1]'R[:,3] +
    λ[6]*R[:,2]'R[:,3]

    dJ = DynamicPolynomials.differentiate(J, [w;g;λ])
    result = solve(dJ)

    reals = realsolutions(result)
    costfun = w-> begin
        sum((A*vec(w[1:9])-[b(w[10:12]) for b in B]).^2)
    end
    # R = tomat(reals[minindex])
    result, costfun
end



result, J = calibForceHomotopy(POSES,forces)
reals = realsolutions(result)
tomat(x) = reshape(x[1:9],3,3) |> copy

Rs = tomat.(reals)
RTR(R) = R'R
minval, minindex = findmin(J.(reals))
R = Rs[minindex]
