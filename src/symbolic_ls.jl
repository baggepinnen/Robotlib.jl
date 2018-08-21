using SymPy
using SymPy: expand, coeff
SymPy.Sym(name, size::Tuple;kwargs...) = Sym[symbols("$name$i$j";kwargs...) for i = 1:size[1], j=1:size[2]]
skewcoords(R) = [R[3,2];R[1,3];R[2,1]]
twistcoords(xi) = [xi[1:3, 4]; skewcoords(xi[1:3, 1:3])]
skew(s) = [0 -s[3] s[2];s[3] 0 -s[1]; -s[2] s[1] 0]
trinv(T) = [T[1:3,1:3].' -T[1:3,1:3].'*T[1:3,4];0 0 0 1]

## Calib NAXP
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





## Calib NAXP with Cayley transform

F = Sym("F",(3,3),real=true)
d = Sym("d",(3,1),real=true)
N =  Sym("N",(3,1),real=true)
PS =  Sym("PS",(3,1),real=true)
PRB =  Sym("PRB",(3,1),real=true)


RA = Sym("RA",(3,3), real=true)
TA = Sym("TA",(3,1), real=true)
g = Sym("g",(3,1), real=true)
u = Sym("u",(3,1), real=true)
G = skew(g)

S = [G u;0 0 0 0]
A = [RA TA;0 0 0 1]
Ps = [PS;1]
Prb = [PRB;1]
N0 = [N;0]

w = [g;u][:]

# s = Ps + inv(Trb)*Pr
# d = -(Ps - inv(Trb)*Pr)
#
# I4 = Matrix{Int}(I,4,4)
# Prb = A*inv(I4-S)*(I4+S)*Ps


eq1 = N0'*[G u;0 0 0 0]*(Ps + trinv(A)*Prb)

eq = eq1[:] # Vector of interesting equations,
# terms of order zeros should be placed on RHS

AA     = Array{Sym}(length(eq), length(w))
BB     = Array{Sym}(length(eq), 1)
AAz    = falses(size(AA))

for i = 1:length(eq)
    for j = 1:length(w)
        AA[i,j] = coeffs(expand(eq[i]),w[j])[2] # Find all coefficients of parameter j in eq i
    end
end

# byhand = (reshape(-repmat((N'A[1:3,1:3]),3,1)',9,1) .*[reshape(repmat(P[1:2],1,3)',6,1);2N])
byhand2 = -skew(N)*[Ps[1:3] (trinv(A)*Prb)[1:3]]
byhand2 = sum(byhand2, dims=2)
byhand2 = [byhand2; 2N]
# @assert all(simplify.(byhand) .== simplify.(byhand2[:]))
@assert all(simplify.(byhand2 .- AA[:]) .== 0)

r = map(AA) do a
    map(N) do n
        coeff(expand(a),n)
    end
end
r = hcat(r...)
