using Robotlib, Test
@testset "Calib Force" begin
N         = 100
Rf        = toOrthoNormal(randn(3,3))
mf        = 1
POSES     = cat([toOrthoNormal(randn(3,3)) for i = 1:N]..., dims=3)
forceoffs = zeros(3)
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*[0, 0, mf*-9.82] - forceoffs) for i = 1:N]...)'

R,m = Robotlib.Calibration.calibForce(POSES,forces; offset=false)

@test R ≈ Rf
@test m ≈ mf


forceoffs = randn(3)
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*[0, 0, mf*-9.82] - forceoffs) for i = 1:N]...)'

R,m,offs = Robotlib.Calibration.calibForce(POSES,forces; offset=true)

@test R ≈ Rf
@test m ≈ mf
@test offs ≈ forceoffs

end
