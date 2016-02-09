using Robotlib
Tworld2base = eye(4)
# """
# `τ = gravity{P}(q::VecOrMat{P}, rm, m, dh::DH, Tworld2base = eye(4))`
# calculates the gravity torques given a joint angles `q`, (distances to center of mass × mass) `rm` and masses `m`. The DH-parameters for the robot must also be provided
# """
function gravity{P}(q::VecOrMat{P}, rm, m, dh::DH, Tworld2base = eye(4))
    Tn  = dh2Tn(dh,q)
    n   = size(Tn,3)-1
    # Find elements that are ≈ 0 and make them zero, smiplifies a lot!
    Tnz = Bool[round(Float64(evalf(subs(Tn[i,j,k],(q[1],1),(q[2],1),(q[3],1),(q[4],1),(q[5],1),(q[6],1),(q[7],1)))),10) == 0 for i=1:4, j=1:4, k=1:n+1]
    Tn[Tnz] = 0
    # Replace the numerical value for π/2 by the SymPi constant PI, this allows for nice simplifications
    Tn = Sym[subs(Tn[i,j,k],("1.5707963267948966",SymPy.PI))  for i=1:4, j=1:4, k=1:n+1]


    Ti = zeros(typeof(q[1]),4,4,n+1)
    Tn = cat(3,Tworld2base,Tn)
    Ti[:,:,1] = Tn[:,:,1]

    for i = 2:n
        Ti[:,:,i] = Ti[:,:,i-1]*Tn[:,:,i]
    end

    gv          = [0, 0, -9.82]
    i           = n
    gi          = Ti[1:3,1:3,i]'*gv
    τ           = Array(Sym,3,n)
    τhat        = Array(Sym,n)
    force       = Array(Sym,3,n)

    τ[:,i]      = skew(rm[:,i]) * gi
    force[:,i]  = gv*m[end]
    τhat[i]     = τ[3,i]

    for i = n-1:-1:1
        Ri      = Tn[1:3,1:3,i+1]
        gi      = Ti[1:3,1:3,i]'*gv # Local frame
        τ[:,i]  = skew(rm[:,i]) * gi + Ri*τ[:,i+1]
        τi      = τ[:,i]
        for k = (i+1):n # This loop formulations seems to be essential, I can not get it to work with accumulated forces
            rk = (trinv(Ti[:,:,i])*Ti[1:4,4,k])[1:3]
            τi += skew(rk)*Ti[1:3,1:3,i]'*force[:,k]
        end

        τhat[i]     = τi[3]
        force[:,i]  = gv*m[i]
    end
    return τhat
end

using Robotlib

using SymPy


# function symbolicGravity(dh::DH)
dh = DHYuMi()
n_joints = size(dh.dhpar,1)

rm  = Sym[symbols("rm$i$j",real=true) for i = 1:3, j=1:n_joints]
m   = Sym[symbols("m$j",real=true) for j=1:n_joints]
q   = Sym[symbols("q$j",real=true) for j=1:n_joints]
tau = Sym[symbols("tau$j",real=true) for j=1:n_joints]


baseAnglesLeft  = [-0.63 , 0.95 , -0.18]
Rbase           = rpy2R(baseAnglesLeft,"xyz")
Tbase           = eye(4)
Tbase[1:3,1:3]  = Rbase

tauhat          = gravity(q, rm, m, dh, Tbase);
w               = [rm[:];m[:]];

AA  = Array(Sym,n_joints, 4n_joints);
AAz = falses(size(AA))

eq  = tauhat;

println("Finding non-zero coefficients")
for i = 1:n_joints
    for j = 1:size(AA,2)
        AA[i,j] = (coeff(expand(eq[i]),w[j])) # This is a nice place to put simplify
        AAz[i,j] = abs(evalf(subs(AA[i,j],(q[1],1),(q[2],1),(q[3],1),(q[4],1),(q[5],1),(q[6],1),(q[7],1)))) < 1e-5
    end
    print(i)
end
println("")


# non_identifiable    = (sum(AA,1) .== 0)[:]
# # non_identifiable = all(AAz,1)[:]
# w                   = w[!non_identifiable];
# AA                  = AA[:,!non_identifiable];
n_joints,n_params   = size(AA)

println("Trying to find expressions common in all terms")
change,res  = cse(AA[:]);
res         = reshape(res,size(AA));

println("Printing results to file")
fid = open("gravityFridaLS.jl","w");
println(fid,"function gravityFridaLS(q)")
for i = 1:n_joints
    println(fid,"q$i = q[$i]")
end
for i = 0:length(change)-1
    sei = change[i+1][2]
    println(fid,"x$i = $sei")
end

println(fid,"A = zeros($n_joints,$n_params)")

for j = 1:n_joints
    for p = 1:n_params
        if res[j,p] != 0 #&& !AAz[j,p]
            println(fid,"A[$j,$p] = $(res[j,p])")
        end
    end
end
println(fid,"return A")
println(fid,"end")
close(fid)
println("Done")
# end
