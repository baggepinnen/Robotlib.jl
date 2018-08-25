N         = 100
Rf        = toOrthoNormal(randn(3,3))
mf        = 1
POSES     = cat([toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
forceoffs = zeros(3)
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*[0, 0, mf*-9.82] - forceoffs) for i = 1:N]...)'

Rf,m = calibForce(POSES,F,m0=0.3; offset=false)
