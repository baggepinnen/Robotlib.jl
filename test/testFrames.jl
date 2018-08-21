using Robotlib, Robotlib.Frames, MAT, Test
R = toOrthoNormal(randn(3,3))
t = randn(3)
f = Frame(R,t)
f = Frame([R t;0 0 0 1])
F = f

path = joinpath(dirname(@__FILE__),"..","src","applications","frames/")
add_frame_name!("SEAM","Weld seam frame")
add_frame_name!("TAB","Table frame")

T_RB_Tm = MAT.matread(path*"T_RB_T.mat")["T_RB_T"]
T_TF_TCPm = MAT.matread(path*"T_TF_TCP.mat")["T_TF_TCP"]
T_T_TABm = MAT.matread(path*"T_T_Table.mat")["T_T_Table"]

T_RB_T = Frame(T_RB_Tm,"RB","T")
T_S_D = Frame(T_TF_TCPm,"S","D")
T_T_TAB = Frame(T_T_TABm,"T","TAB")

cloud_seam = readcloud(path*"CloudSeam_edge.txt")
plane_seam = readplane(path*"PlaneSeam_edge.txt")
cloud_seam_projected = project(plane_seam,cloud_seam)
line_seam = fitline(cloud_seam_projected)

T_T_SEAM = framefromfeatures(("z+",line_seam),("y-",plane_seam),cloud_seam_projected[1],"SEAM")
T_RB_SEAM = T_RB_T*T_T_SEAM
T_RB_TAB = T_RB_T*T_T_TAB
T_TAB_SEAM = inv(T_T_TAB)*T_T_SEAM

cloud_seam_RB = T_RB_T*cloud_seam
cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
plane_seam_RB = T_RB_T*plane_seam
line_seam_RB = T_RB_T*line_seam



# plot(Frame(I4,"RB","U"),200)
# plot!(cloud_seam_RB,label="Cloud seam_RB")
# plot!(cloud_seam_projected_RB,label="Cloud seam_RB projected")
# plot!(line_seam_RB,500,label="Line seam")
# plot!(plane_seam_RB,200,label="Plane seam")
# plot!(T_RB_SEAM,200)
# plot!(T_RB_TAB,200)
# gui()


T = convert(Matrix, f)
@test convert(Array,f) == T
@test promote(T,f) == (T,T)
@test promote_type(Matrix{Float64}, Frame{Float64}) == Matrix{Float64}
