"""
`calibForce(POSES, F, m0=0.3; offset=true)`
If result is bad, check if you send data in correct form;
POSES is always the position of the tool frame from the robot FK
cell array with 4x4 transformation matrices
F is a Nx3 vector of forces (accepts Nx6 matrix with torques also)
usage `Rf*force[i,1:3]' + forceoffs = POSES[1:3,1:3,i]'*[0, 0, mf*-9.82]`
Bagge
"""
function calibForce(POSES,F,m0=0.3; offset=true, useCVX = false)

    N = size(POSES,3)

    forceoffs = zeros(3)
    g = 9.82
    mg = m0*g # Initial guess is m0, the method accepts 5 orders of magnitude error at least

    I = I3
    B = Array{Float64}(3N)
    B2 = Array{Float64}(3N)
    if offset # Include offset
        A = Array{Float64}(3N,12)
        for i = 1:N
            RA = POSES[1:3,1:3,i]'
            At = [F[i,1]*I   F[i,2]*I    F[i,3]*I   I]
            b =  -RA[1:3,3]
            A[3(i-1)+1:3i,:] = At
            B[3(i-1)+1:3i] = b
        end
        info("Condition number of Gram matrix: ", cond(A'A))
        w = A\B.*mg
        if useCVX
            w = Variable(12)
            problem = minimize(norm(A*w-B.*mg,1))
            solve!(problem)
        end
        Rf = reshape(w[1:9],3,3)
        info("Determinant of Rf before orthonormalization: ", det(Rf)," (Should be close to 1, but don't be afraid if it's not, should however be positive!)")
        det(Rf) < 0 && @warn("det(Rf) < 0, left handed coordinate system?")
        toOrthoNormal!(Rf)
        for ii = 1:2
            #        forceoffs = w(10:12)
            A2 = Array{Float64}(3N,4)
            for i = 1:N
                RA = POSES[1:3,1:3,i]'
                At = [RA[:,3] I]
                b =  -Rf*F[i,1:3]
                A2[3(i-1)+1:3i,:] = At
                B2[3(i-1)+1:3i] = b
            end
            w = A2\B2
            if useCVX
                w = Variable(4)
                problem = minimize(norm(A2*w-B2,1))
                solve!(problem)
            end
            mg = w[1]
            forceoffs = w[2:4]
        end
    else # Do not include offset
        A = Array{Float64}(3N,9)
        for i = 1:N
            RA = POSES[1:3,1:3,i]'
            At = [F[i,1]*I F[i,2]*I F[i,3]*I]
            b =  -RA[1:3,3]
            A[3(i-1)+1:3i,:] = At
            B[3(i-1)+1:3i] = b
        end
        info("Condition number of Gram matrix: ", cond(A'A))
        w = A\B.*mg
        for ii = 1:4
            Rf = reshape(w[1:9],3,3)
            info("Determinant of Rf before orthonormalization: ", det(Rf))
            det(Rf) < 0 && @warn("det(Rf) < 0, left handed coordinate system?")
            toOrthoNormal!(Rf)
            A2 = Array{Float64}(3N,1)
            B2 = Array{Float64}(3N)
            for i = 1:N
                RA = POSES[1:3,1:3,i]'
                At = RA[:,3]
                b =  -Rf*F[i,1:3]
                A2[3(i-1)+1:3i,:] = At
                B2[3(i-1)+1:3i] = b
            end
            mg = (A2\B2)[]

        end
    end
    m = mg/g
    if m < 0
        @warn("Estimated mass is negative, left handed coordinate system for sensor?")
    end
    if offset
        return Rf,m,forceoffs
    else
        return Rf,m
    end
end



# function calibForceCVX(T,f,m=0.3)
#
#     Rfv = Variable(3,3)
#     offsetv = Variable(3)
#     J = 0.0
#     m = 0.2205
#     for i = 1:10:N
#         J += sum(abs(Rfv*f[i,1:3]' + offsetv - T[1:3,1:3,i]'*[0, 0, m*-9.82]))
#     end
#     problem = minimize(J)
#     solve!(problem, SCSSolver(verbose=true,max_iters=10000))
#     Rf = Rfv.value
#     offset = Rf\offsetv.value
#     Rf = toOrthoNormal(Rf)
#     mv = Variable(1)
#     J = 0.0
#     for i = 1:10:N
#         J += sum(abs(Rf*(f[i,1:3]' + offset) - T[1:3,1:3,i]'*[0, 0, mv*-9.82]))
#     end
#     problem = minimize(J)
#     solve!(problem, SCSSolver(verbose=true,max_iters=10000))
#     m = mv.value
#
#     model = JuMP.Model(solver=IpoptSolver())
#     # model = JuMP.Model(solver=IpoptSolver(constr_viol_tol = 10., dual_inf_tol = 10., compl_inf_tol=10.))
#     JuMP.@defVar(model, g[1:3])
#     JuMP.@defVar(model, offset2[1:3])
#     JuMP.@defVar(model, Rf2[1:3,1:3])
#     # JuMP.setValue(Rf2,Rf)
#     JuMP.setValue(g,[0, 0, m*-9.82])
#     JuMP.setValue(offset2,offsetv.value)
#     J = 0.0
#
#     function absv(v::JuMP.Variable)
#         JuMP.@defVar(v.m, aux >= 0)
#         JuMP.@addConstraint(v.m, aux >= v)
#         JuMP.@addConstraint(v.m, aux >= -v)
#         return aux
#     end
#
#     function sumabsv(m,v::VecOrMat)
#         N = length(v)
#         JuMP.@defVar(m, aux[1:N] >= 0)
#         for i = 1:N
#             JuMP.@addConstraint(m, aux[i] >= v[i])
#             JuMP.@addConstraint(m, aux[i] >= -v[i])
#         end
#         return sum(aux)
#     end
#
#     for i = 1:10:N
#         X = Rf2*f[i,1:3]' + offset - T[1:3,1:3,i]'*g
#         J += sumabsv(model,X)
#     end
#
#     for i = 1:3
#         for j = 1:3
#             if i == j
#                 JuMP.@addConstraint(model, Rf2[1:3,i]'*Rf2[1:3,i] .== 1)
#                 # JuMP.@addNLConstraint(model, Rf2'Rf2 == I)
#                 # JuMP.@addNLConstraint(model, sum{Rf2[i,k]*Rf2[k,j], k=1:3} == 1)
#             else
#                 JuMP.@addConstraint(model, Rf2[1:3,i]'*Rf2[1:3,j] .== 0)
#                 # JuMP.@addNLConstraint(model, sum{Rf2[i,k]*Rf2[k,j], k=1:3} == 0)
#             end
#
#         end
#     end
#     JuMP.@setObjective(model, Min, J)
#     JuMP.solve(model)
#     Rf2 = JuMP.getValue(Rf2)
#     g = JuMP.getValue(g)
#
#     err = cat(2,[Rf*(f[i,1:3]' + offset) - T[1:3,1:3,i]'*[0, 0, m*-9.82] for i = 1:N]...)'
#     err2 = cat(2,[Rf2*(f[i,1:3]' + offset) - T[1:3,1:3,i]'*g for i = 1:N]...)'
#     plot(f[:,1:3],lab="Force")
#     plot!(err,l=:dash,lab="Error")
#     plot!(err,l=:dash,lab="Error2")
#     println("Error:  ", round(rms(err),4))
#     println("Error2: ", round(rms(err2),4))
#
# end
