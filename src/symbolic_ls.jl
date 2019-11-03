# This script contains some code for symbolic manipulation of linear least-squares problems. It can be used to derive or verify matrix expressions associated with problems on the form Ax = y
using SymPy
using SymPy: expand, coeff
SymPy.Sym(name, size::Tuple;kwargs...) = Sym[symbols("$name$i$j";kwargs...) for i = 1:size[1], j=1:size[2]]
skewcoords(R) = [R[3,2];R[1,3];R[2,1]]
twistcoords(xi) = [xi[1:3, 4]; skewcoords(xi[1:3, 1:3])]
skew(s) = [0 -s[3] s[2];s[3] 0 -s[1]; -s[2] s[1] 0]
trinv(T) = [T[1:3,1:3]' -T[1:3,1:3]'*T[1:3,4];0 0 0 1]

## Calib NAXP ==============================================================================
# N' * H*Pf = norm(N)^2


F = Sym("F",(3,3),real=true)
d = Sym("d",(3,1),real=true)
N =  Sym("N",(3,1),real=true)
P =  Sym("P_",(3,1),real=true)


RA = Sym("RA",(3,3), real=true)
TA = Sym("TA",(3,1), real=true)
RX = Sym("RX",(3,3), real=true)
TX = Sym("TX",(3,1), real=true)
X = [RX TX;0 0 0 1]
A = [RA TA;0 0 0 1]
Ps = [P;1]
N0 = [N;0]

xh = symbols("xh", real=true)
yh = symbols("yh", real=true)


# H = F'[[1 0;0 1;0 0], -d];

w = X[1:3,[1, 2, 4]]
w = w[:]
Prb = A*(X*Ps)
eq1 = N'*Prb[1:3]
eq = eq1[:] # Vector of interesting equations,
# terms of order zeros should be placed on RHS

AA     = Array{Sym}(length(eq), length(w))
BB     = Array{Sym}(length(eq), 1)
AAz    = falses(size(AA))

for i = 1:length(eq)
    for j = 1:length(w)
        AA[i,j] = coeff(expand(eq[i]),w[j]) # Find all coefficients of parameter j in eq i
    end
end

println("Trying to find expressions common in all terms")
change,res  = cse(AA[:]);
res         = reshape(res,size(AA));


byhand = (reshape(repmat((N'A[1:3,1:3]),3,1)',9,1) .*[reshape(repmat(P[1:2],1,3)',6,1);1;1;1])
byhand2 = [N'A[1:3,1:3]*Ps[1] N'A[1:3,1:3]*Ps[2] N'A[1:3,1:3]]


@assert all(simplify.(byhand) .== simplify.(byhand2[:]))
@assert all(simplify.(byhand) .== simplify.(AA[:]))
