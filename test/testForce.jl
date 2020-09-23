using Robotlib, Test, Random
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

R = Robotlib.Calibration.calibForce2(POSES,forces,mf)
@test R ≈ Rf

forceoffs = randn(3)
forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*[0, 0, mf*-9.82] - forceoffs) for i = 1:N]...)'

R,m,offs = Robotlib.Calibration.calibForce(POSES,forces; offset=true)

@test R ≈ Rf
@test m ≈ mf
@test offs ≈ forceoffs


end



# using Plots, DSP, Robotlib, Random, LaTeXStrings
# pgfplots()
# σr = exp10.(range(-4,stop=1,length=1000))
# res = map(σr) do σ
#
#     N         = 100
#     Rf        = expω(1randn(3))
#     mf        = 1
#     POSES     = cat([expω(randn(3)) for i = 1:N]..., dims=3)
#     forceoffs = zeros(3)
#     forces    = hcat([Rf'*(POSES[1:3,1:3,i]'*[0, 0, mf*-9.82] - forceoffs) for i = 1:N]...)'
#     forces  .+= σ*randn(size(forces))
#
#     R1,m = Robotlib.Calibration.calibForce(POSES,forces,1.1mf; offset=false)
#
#     R2 = Robotlib.Calibration.calibForce2(POSES,forces,1.1mf)
#     sum(abs2,R1-Rf), sum(abs2,R2-Rf)
# end
#
# e1,e2 = getindex.(res,1),getindex.(res,2)
# plot(σr,filtfilt(ones(20),[20],[e1 e2]), xscale=:log10, yscale=:log10, lab=["Relaxation" "Cayley"], xlabel=L"$\sigma$ [N]", ylabel=L"$||R_e||_F^2$",legend=:bottomright)
# hline!([3, 3], lab="", l=:dash)
# savefig2("/home/fredrikb/phdthesis/calibpaper/figs/calibforce2.tex")
