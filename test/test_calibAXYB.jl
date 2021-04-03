using Robotlib
using Robotlib: T2t


T_RB_T0, T_TF_S0, POSES, MEASURED, LEDs = Robotlib.simulateCalibration_AXYB(100, σt = 0.1)

T_RB_T, T_TF_S, normDist = calibAXYB(POSES, MEASURED, refine=false)
@test T_RB_T ≈ T_RB_T0 rtol = 1e-6
@test T_TF_S ≈ T_TF_S0 rtol = 1e-6
@test length(LEDs) == 3
@test size(LEDs[1]) == (3, 100)
L = leds_to_tf(LEDs)
@test LEDs[1][:,1] ≈ T2t(MEASURED[1])
@test T2t(L[1]) ≈ T2t(MEASURED[1])
@test maximum(norm.(L .- MEASURED)) < 1e-6

## With noise

T_RB_T0, T_TF_S0, POSES, MEASURED, LEDs = Robotlib.simulateCalibration_AXYB(100, σt = 0.1, σy = 0.05)

L = leds_to_tf(LEDs)
T_RB_T, T_TF_S, normDist = calibAXYB(POSES, L, refine=false)
@test T2t(T_RB_T) ≈ T2t(T_RB_T0) rtol = 3e-1
@test T2t(T_TF_S) ≈ T2t(T_TF_S0) rtol = 3e-1

@test Rangle(T_RB_T, T_RB_T0) < deg2rad(5)
@test Rangle(T_TF_S, T_TF_S0) < deg2rad(25)

T_RB_T2, T_TF_S2, normDist2 = calibAXYB(POSES, L, refine=true)
@test T2t(T_RB_T2) ≈ T2t(T_RB_T0) rtol = 3e-1
@test T2t(T_TF_S2) ≈ T2t(T_TF_S0) rtol = 3e-1

@test Rangle(T_RB_T2, T_RB_T0) < deg2rad(5)
@test Rangle(T_TF_S2, T_TF_S0) < deg2rad(25)

# @test norm(T2t(T_RB_T2) - T2t(T_RB_T0)) < norm(T2t(T_RB_T) - T2t(T_RB_T0)) # these are not stable
# @test norm(T2t(T_TF_S2) - T2t(T_TF_S0)) < norm(T2t(T_TF_S) - T2t(T_TF_S0)) # these are not stable

# @test Rangle(T_RB_T2, T_RB_T0) < Rangle(T_RB_T, T_RB_T0) # these are not stable
# @test Rangle(T_TF_S2, T_TF_S0) < Rangle(T_TF_S, T_TF_S0) # these are not stable

@test sum(normDist2) < sum(normDist)