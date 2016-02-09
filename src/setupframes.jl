using Frames
import MAT
function setupframes(path)
	#include("/home/fredrikb/robotlib/frames/Frames.jl")
	# path = "/home/fredrikb/work/flexifab/frames/"


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

	MAT.matwrite(path*"T_TAB_SEAM.mat",["T_TAB_SEAM" => T_TAB_SEAM.T])
	MAT.matwrite(path*"T_T_SEAM.mat",["T_T_SEAM" => T_T_SEAM.T])
	MAT.matwrite(path*"T_RB_TAB.mat",["T_RB_TAB" => T_RB_TAB.T])
	println("Wrote T_TAB_SEAM, T_T_SEAM, T_RB_TAB to files in $path")

	cloud_seam_RB = T_RB_T*cloud_seam
	cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
	plane_seam_RB = T_RB_T*plane_seam
	line_seam_RB = T_RB_T*line_seam



	plotframe(Frame(eye(4),"RB","U"),200, label=true)

	plotpoints(cloud_seam_RB,"bx")
	plotpoints(cloud_seam_projected_RB,"rx")
	plotline(line_seam_RB,"r",500,label="Line seam")
	plotplane(plane_seam_RB,"b",200,label="Plane seam")
	plotframe(T_RB_SEAM,200, label=true)
	plotframe(T_RB_TAB,200, label=true)

	PyPlot.xlabel("x")
	PyPlot.ylabel("y")
	PyPlot.zlabel("z")
	PyPlot.axis("scaled")


	if true
		add_frame_name!("SEAM2","Frame for second part of Weld seam")
		cloud_seam2 = readcloud(path*"CloudSeam2_edge.txt")
		cloud_seam2_projected = project(plane_seam,cloud_seam2)
		line_seam2 = fitline(cloud_seam2_projected)

		T_T_SEAM2 = framefromfeatures(("z+",line_seam2),("y-",plane_seam),cloud_seam2_projected[1],"SEAM2")
		T_RB_SEAM2 = T_RB_T*T_T_SEAM2
		T_TAB_SEAM2 = inv(T_T_TAB)*T_T_SEAM2;
		MAT.matwrite("T_TAB_SEAM2.mat",["T_TAB_SEAM2" => T_TAB_SEAM2.T])
		MAT.matwrite("path*"T_T_SEAM2.mat",["T_T_SEAM2" => T_T_SEAM2.T])
		println("Wrote T_TAB_SEAM2, T_T_SEAM2 to files in $path")


		cloud_seam2_RB = T_RB_T*cloud_seam2
		cloud_seam2_projected_RB = T_RB_T*cloud_seam2_projected
		plane_seam2_RB = T_RB_T*plane_seam
		line_seam2_RB = T_RB_T*line_seam2

		plotpoints(cloud_seam2_RB,"bx")
		plotpoints(cloud_seam2_projected_RB,"rx")
		plotline(line_seam2_RB,"r",500,label="Line seam2")
		plotplane(plane_seam2_RB,"b",200,label="Plane seam2")
		plotframe(T_RB_SEAM2,200, label=true)
	end

end
