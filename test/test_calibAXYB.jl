using Robotlib

T_RB_T0, T_TF_S0, POSES, MEASURED, LEDs = Robotlib.simulateCalibration_AXYB(100, σt = 0.1)

T_RB_T, T_TF_S, normDist = calibAXYB(POSES, MEASURED)

@test T_RB_T ≈ T_RB_T0 rtol = 1e-6
@test T_TF_S ≈ T_TF_S0 rtol = 1e-6

@test length(LEDs) == 3
@test size(LEDs[1]) == (3, 100)

L = leds_to_tf(LEDs)

@test LEDs[1][:,1] ≈ T2t(MEASURED[1])
@test T2t(L[1]) ≈ T2t(MEASURED[1])


@test maximum(norm.(L .- MEASURED)) < 1e-6