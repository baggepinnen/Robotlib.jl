"""
Frames is a module to work with coordinate systems. The central type is the `Frame` with fields `R` and `t` for the rotation matrix and translation vector respectively. Other inportnat important types include `Point, Points, Plane` and `Line`

All geometric object types have a reference frame associated with them. When a geometric object is created, the reference frame is specified, e.g., `T_RB_T = Frame(T_RB_Tm,"RB","T")` is a transformation between frame `RB` and frame `T`.
When transofrmations are made, the reference frames are updated automatically, e.g., `T_RB_SEAM = T_RB_T*T_T_SEAM`

Usage:
```julia
using Robotlib.Frames
path = Pkg.dir("Robotlib","src","applications","frames/")
add_frame_name!("SEAM","Weld seam frame")
add_frame_name!("TAB","Table frame")

T_RB_Tm    = MAT.matread(path*"T_RB_T.mat")["T_RB_T"]
T_TF_TCPm  = MAT.matread(path*"T_TF_TCP.mat")["T_TF_TCP"]
T_T_TABm   = MAT.matread(path*"T_T_Table.mat")["T_T_Table"]

T_RB_T     = Frame(T_RB_Tm,"RB","T")
T_S_D      = Frame(T_TF_TCPm,"S","D")
T_T_TAB    = Frame(T_T_TABm,"T","TAB")

cloud_seam = readcloud(path*"CloudSeam_edge.txt")
plane_seam = readplane(path*"PlaneSeam_edge.txt")
cloud_seam_projected = project(plane_seam,cloud_seam)
line_seam  = fitline(cloud_seam_projected)

T_T_SEAM      = framefromfeatures(("z+",line_seam),("y-",plane_seam),cloud_seam_projected[1],"SEAM")
T_RB_SEAM     = T_RB_T*T_T_SEAM
T_RB_TAB      = T_RB_T*T_T_TAB
T_TAB_SEAM    = inv(T_T_TAB)*T_T_SEAM

cloud_seam_RB = T_RB_T*cloud_seam
cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
plane_seam_RB = T_RB_T*plane_seam
line_seam_RB  = T_RB_T*line_seam

plot(Frame(eye(4),"RB","U"),200, label=true)
plot!(cloud_seam_RB)
plot!(cloud_seam_projected_RB)
plot!(line_seam_RB,500,label="Line seam")
plot!(plane_seam_RB,200,label="Plane seam")
plot!(T_RB_SEAM,200, label=true)
plot!(T_RB_TAB,200, label=true)

xlabel!("x")
ylabel!("y")
zlabel!("z")

MAT.matwrite(path*"T_TAB_SEAM.mat",["T_TAB_SEAM" => T_TAB_SEAM.T])
MAT.matwrite(path*"T_T_SEAM.mat",["T_T_SEAM" => T_T_SEAM.T])
MAT.matwrite(path*"T_RB_TAB.mat",["T_RB_TAB" => T_RB_TAB.T])
println("Wrote T_TAB_SEAM, T_T_SEAM, T_RB_TAB to files in \$path")
```
"""
module Frames

using Plots, FixedSizeArrays
default(markersize=1)
export Frame, Point, Plane, Points, Line, GeometricObject, add_frame_name!
export readcloud, readTmatrix, readplane, fitline, fitplane, framefromfeatures, project
export plot3Dsmart, display, show, print
export inv, *,+,-,/,\,transpose,ctranspose, dot

import Base: det, print, zeros, length, size, getindex, setindex!, convert, promote_rule, push!, show, display, start, next, done, +, *, ⋅, .*, /, ./, -, ×, transpose, ctranspose, \, inv, dot
# using LaTeXStrings
import Robotlib: T2R, T2t

typealias Vect{Ty} Union{AbstractVector{Ty}, Vec{3,Ty}}
typealias Matr{Ty} Union{AbstractMatrix{Ty}, Mat{3,3,Ty}}

# Type definitions --------------------------------------------
# -------------------------------------------------------------
abstract GeometricObject

<<<<<<< HEAD
type Point <: GeometricObject
    p::Vector{Float64}
    A::String
    function Point(p::Vector{Float64} = zeros(3), A::String="U")
        checkframe(A,"U")
        new(p,A)
    end
end

type Points <: GeometricObject
    p::Vector{Point}
    A::String
    function Points(p::Vector{Point}, A::String="U")
        checkframe(A,"U")
        new(p,A)
    end

    function Points(p::Array{Float64,2}, A::String="U")
        checkframe(A,"U")
        N = size(p,1)
        points = [Point(squeeze(p[i,:]',2),A) for i in 1:N]
        new(points,A)
    end
    function Points(A::String="U")
        new(Point[],A)
    end
=======
type Point{Ty} <: GeometricObject
    p::Vec{3,Ty}
    A::String
end
function Point{Ty}(p::Vect{Ty} = zeros(3), A::String="U")
    checkframe(A,"U")
    Point{Ty}(p,A)
end

type Points{Ty} <: GeometricObject
    p::Vector{Point{Ty}}
    A::String
end
function Points{Ty}(p::AbstractVector{Point{Ty}}, A::String="U")
    checkframe(A,"U")
    Points{Ty}(p,A)
end

function Points{Ty}(p::AbstractMatrix{Ty}, A::String="U")
    checkframe(A,"U")
    N = size(p,1)
    points = [Point{Ty}(squeeze(p[i,:]',2),A) for i in 1:N]
    Points{Ty}(points,A)
end
function Points(A::String="U")
    Points{Float64}(Point{Float64}[],A)
>>>>>>> 4e7648502043a3c813a37584688307a109a586bd
end

push!(points::Points, p::Point) = push!(points.p,p)

getindex(points::Points,i) = points.p[i]
function setindex!(points::Points, p::Point, i)
    points.p[i] = p
end
getindex(point::Point,i) = point.p[i]
function setindex!(point::Point, p::Float64, i)
    point.p[i] = p
end

size(points::Points) = size(points.p)[1]

length(points::Points) = length(points.p)

Base.start(points::Points) = Base.start(points.p)
Base.next(points::Points, state) = Base.next(points.p,state)
Base.done(points::Points, state) = Base.done(points.p,state)

function convert(::Type{Matrix{Float64}}, p::Points)
    N = length(p)
    points = zeros(N,3)
    for i = 1:N
        points[i,:] = Vector(p[i].p)
    end
    return points
end

function zeros(points::Points)
    p = Point()
    N = length(points)
    newpoints = [p for i in 1:N]
    Points(newpoints,points.A)
end

type Frame{Ty} <: GeometricObject
    R::Mat{3,3,Ty}
    t::Vec{3,Ty}
    A::String
    B::String
end
function Frame{Ty}(T::Matr{Ty}, A::String="U",B::String="U")
    checkframe(A,B)
    f = Frame{Ty}(T2R(T),T2t(T),A,B)
    add_frame!(f)
    f
end
function Frame{Ty}(R::Matr{Ty}, t::AbstractArray{Ty}, A::String="U", B::String="U")
    checkframe(A,B)
    f = Frame{Ty}(R,t,A,B)
    add_frame!(f)
    f
end
Frame() = Frame{Float64}(eye(3),zeros(3),"U","U")

<<<<<<< HEAD
    T::Matrix{Float64}
    A::String
    B::String
    function Frame(T::Matrix, A::String="U",B::String="U")
        checkframe(A,B)
        f = new(T,A,B)
        add_frame!(f)
        f
    end
    function Frame(R::Matrix, t::Array, A::String="U", B::String="U")
        checkframe(A,B)
        f = new([R t; 0 0 0 1],A,B)
        add_frame!(f)
        f
    end
    function Frame(); f = new(eye(4),"U","U"); add_frame!(f);  f end
=======
# Frame(R::Matr, t::AbstractArray, A::String="U", B::String="U") = Frame{eltype(R)}(R, t, A, B)
>>>>>>> 4e7648502043a3c813a37584688307a109a586bd


ref_frames = Set{Frame}()
frame_map = Dict{String,String}()
frame_names = Dict("U" => "Undefined", "RB"=>"Robot Base","TF"=>"Tool Flange","TCP"=>"Tool Center Point","T"=>"Tracker","D"=>"Dynamic diods","S"=>"Sensor")
add_frame_name!(key::String,value::String) = setindex!(frame_names,value,key)
frame_name_exists(ref::String) = ref ∈ keys(frame_names)
function add_frame!(f::Frame)
    push!(ref_frames,f)
    setindex!(frame_map,f.B,f.A)
end

function checkframe(A,B)
    if !frame_name_exists(A)
        error("Frame \"$A\" does not exist, add to dictionary using add_ref_frame()")
    elseif !frame_name_exists(B)
        error("Frame \"$B\" does not exist, add to dictionary using add_ref_frame()")
    end
end

checkframes(x::GeometricObject, y::GeometricObject) = x.A != y.A &&  error("The two geometric objects do not have the same reference frame ($x -> $(x.A), $y -> $(y.A))")

function display(f::Frame)
    println(round(F2T(f),4))
    show(string(frame_names[f.A])*"->"*string(frame_names[f.B]))
end
print(f::Frame) = "\$T_{$(f.A)}^{$(f.B)}\$"
show(f::Frame) = display(f)

<<<<<<< HEAD
type Plane <: GeometricObject
    n::Point
    r::Point
    A::String
    Plane(n,r,A="U") = new(n,r,A)
    Plane(n::Vector,r::Vector,A="U") = new(Point(n,A),Point(r,A),A)
    Plane(points::Points) = fitplane(points)
=======
type Plane{Ty} <: GeometricObject
    n::Point{Ty}
    r::Point{Ty}
    A::String
>>>>>>> 4e7648502043a3c813a37584688307a109a586bd
end
Plane{Ty}(n::Vect{Ty},r::Vect{Ty},A="U") = Plane{Ty}(n,r,A)
Plane{Ty}(n::Vect{Ty},r::Vect{Ty},A="U") = Plane{Ty}(Point{Ty}(n,A),Point{Ty}(r,A),A)
Plane(points::Points) = fitplane(points)

<<<<<<< HEAD
type Line <: GeometricObject
    v::Point
    r::Point
    A::String
    Line(v,r,A="U") = new(v,r,A)
    Line(points::Points) = fitline(points)
=======
type Line{Ty} <: GeometricObject
    v::Point{Ty}
    r::Point{Ty}
    A::String
>>>>>>> 4e7648502043a3c813a37584688307a109a586bd
end
Line{Ty}(v::Vect{Ty},r::Vect{Ty},A="U") = Line{Ty}(v,r,A)
Line(points::Points) = fitline(points)




# Plot functions ----------------------------------------------
# -------------------------------------------------------------

mat2tup(m) = m[:,1][:], m[:,2][:], m[:,3][:]

@recipe function plotframe(f::Frame, length=1.0)#; annotation=false)
    #Plots a frame using XYZ-RGB
    o = F2t(f) |> Vector
    x = Rx(f) |> Vector
    y = Ry(f) |> Vector
    z = Rz(f) |> Vector
    seriestype := :path3d
    label     := ""
    @series begin
        data = [o o+x*length]'
        seriescolor := :red
        data |> mat2tup
    end
    @series begin
        data = [o o+y*length]'
        seriescolor := :green
        data |> mat2tup
    end
    # if annotation && f.A != ""
    #     warn("Annotations not currently supported for 3d plots")
    #     # po = o-length/4*(x+y+z)
    #     # annotations := (po[1],po[2],po[3],print(f))
    # end
    # delete!(d,:annotation)
    @series begin
        data = [o o+z*length]'
        seriescolor := :blue
        data |> mat2tup
    end
end



@recipe function plotline(l::Line, length=1.0)#; annotation="")
    coords = [Vector((l.r-length*l.v).p)  Vector((l.r+length*l.v).p)]'
    @series begin
        seriestype := :scatter3d
        ([l.r[1]],[l.r[2]],[l.r[3]])
    end
    # if annotation != ""
    #     warn("Annotations not currently supported for 3d plots")
    #     # po = l.r.p-length/4*normalized(l.r.p)
    #     # annotations := (po[1],po[2],po[3],d[:annotation])
    # end
    # delete!(d,:annotation)
    @series begin
        seriestype := :path3d
        coords |> mat2tup
    end
end


@recipe function plotplane(p::Plane, length=1.0)
    coords = [Vector((p.r).p) Vector((p.r+length*p.n).p)]'
    @series begin
        seriestype := :scatter3d
        markershape --> :utriangle
        ([p.r[1]],[p.r[2]],[p.r[3]])
    end
    @series begin
        seriestype := :path3d
        coords |> mat2tup
    end
    # if label != ""
    #     po = p.r.p-length/4*p.n.p
    #     annotate!(po[1],po[2],po[3],label)
    # end
end

@recipe function plotpoints(p::Points)
    seriestype := :scatter3d
    convert(Matrix{Float64},p) |> mat2tup
end



# Frame functions  --------------------------------------------
# -------------------------------------------------------------

*(p::Point, c::Real) = Point(c.*p.p, p.A)
*(c::Real, p::Point) = *(p, c)
.*(p::Point, c::Real) = *(p, c)
.*(c::Real, p::Point) = *(p, c)
*(p::Point, R::Matr) = error("Post multiplication of ::Point with ::Matrix not supported")
*(R::Matr, p::Point) = Point(R*p.p, p.A)
/(p::Point, c::Real) = Point(p.p./c, p.A)
/(c::Real,p::Point) = /(p,c)
./(p::Point, c::Real) = Point(p,c)
./(c::Real,p::Point) = /(p,c)
+(p1::Point, p2::Point) = Point(p1.p + p2.p, p1.A)
-(p1::Point, p2::Point) = Point(p1.p - p2.p, p1.A)
+(p1::Vect, p2::Point) = p1 + p2.p
-(p1::Vect, p2::Point) = p1 - p2.p
+(p2::Point, p1::Vect) = Point(p1 + p2.p, p2.A)
-(p2::Point, p1::Vect) = Point(p1 - p2.p, p2.A)
⋅(p1::Point, p2::Point) = p1.p ⋅ p2.p
×(p1::Point, p2::Point) = Point(p1.p × p2.p, p1.A)

dot{Ty}(a::Vec{3,Ty}, b::AbstractVector{Ty}) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
dot{Ty}(a::AbstractVector{Ty}, b::Vec{3,Ty}) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

transpose(p::Point) = p.p'
ctranspose(p::Point) = p.p'
convert(::Type{Vector}, p::Point) = p.p
convert{T1,T2,n}(::Type{Array{T1,n}}, f::Frame{T2}) = [Matrix{T1}(f.R) Vector{T1}(f.t); 0 0 0 1]
convert{T2}(::Type{Matrix}, f::Frame{T2}) = convert(Matrix{T2},f)
convert{T2}(::Type{Array}, f::Frame{T2}) = convert(Matrix{T2},f)
promote_rule{T1,T2,n}(::Type{Array{T1,n}}, f::Type{Frame{T2}} ) = Array{promote_type(T1,T2),2}

*(x::Frame, y::Frame) = Frame(x.R*y.R, x.t + x.R*y.t , x.A, y.B)
function *(f::Frame, p::Point)
    f.B != p.A &&  error("Reference frames does not match between $f and $p ($(f.B) != $(p.A))")
    pv = f.R*p.p + f.t
    Point(pv,f.A)
end
function *(f::Frame, p::Plane)
    f.B != p.A && error("Reference frames does not match between $f and $p ($(f.B) != $(p.A))")
    n = F2R(f)*p.n
    r = f*p.r
    λ = (r⋅n)/(n⋅n)
    r = λ*n
    Plane(n,r,f.A)
end
function *(f::Frame, l::Line)
    f.B != l.A && error("Reference frames does not match between $f and $l ($(f.B) != $(l.A))")
    v = F2R(f)*l.v
    r2 = f*l.r
    r3 = v×r2
    r = r3×v
    Line(v,r,f.A)
end
function *(f::Frame, points::Points)
    f.B != points.A &&  error("Reference frames does not match between $f and $p ($(f.B) != $(p.A))")
    newpoints = deepcopy(points)
    newpoints.A = f.A
    for i in 1:length(newpoints)
        newpoints[i] = f*newpoints[i]
    end
    newpoints
end
*(p::Point, f::Frame) = error("Post multiplication of ::Point with ::Frame not supported")
tinv(T)       = [T[1:3,1:3]' -T[1:3,1:3]'*T[1:3,4]; 0 0 0 1]
<<<<<<< HEAD
inv(x::Frame) = Frame(tinv(x.T),x.B,x.A)
\(x::Frame,y::Frame) = Frame(tinv(x.T)*y.T,x.B,y.B)
/(x::Frame,y::Frame) = Frame(x.T*tinv(y.T),x.A,y.A)
F2t(F::Frame) = F.T[1:3,4]
F2R(F::Frame) = F.T[1:3,1:3]
Rx(F::Frame)  = F.T[1:3,1]
Ry(F::Frame)  = F.T[1:3,2]
Rz(F::Frame)  = F.T[1:3,3]
det(F::Frame) = det(F.T[1:3,1:3])
normalized(v::Vector{Float64}) = v/norm(v)
function normalize!(v::Vector{Float64}) v/=norm(v); end
=======
function inv{T}(x::Frame{T}) # TODO: can this implementation be faster?
    R = x.R'
    Frame{T}(R, -R*x.t, x.B,x.A)
end
\{T}(x::Frame{T},y::Frame{T}) = Frame{T}(inv(x)*y,x.B,y.B)
/{T}(x::Frame{T},y::Frame{T}) = Frame{T}(x*inv(y),x.A,y.A)
F2t(F::Frame) = F.t
F2R(F::Frame) = F.R
F2T(F::Frame) = [Array(F.R) Array(F.t); 0 0 0 1]
Rx(F::Frame)  = F.R[1:3,1] |> Vec
Ry(F::Frame)  = F.R[1:3,2] |> Vec
Rz(F::Frame)  = F.R[1:3,3] |> Vec
det(F::Frame) = det(F.R)
normalized(v::Vect) = v/norm(v)
function normalize!(v::Vect) v/=norm(v); end
>>>>>>> 4e7648502043a3c813a37584688307a109a586bd
normalized(p::Point) = Point(p.p/norm(p.p), p.A)
function normalize!(p::Point) p.p/=norm(p.p); end

function ⊥(R::Matr)
    (U,S,V) = svd(R)
    a = diagm([1, 1, sign(det(U*V'))])
    return U*a*V'
end

⊥(F::Frame) = Frame(⊥(F2R(F2)),F.t,F.A,F.B)


# Projections and fitting  ------------------------------------
# -------------------------------------------------------------

function center{Ty}(points::Points{Ty})
    N = size(points)
    cog = zeros(Ty,3)
    for i = 1:N
        cog += points[i]
    end
    cog /= N

    cpoints = zeros(Ty,length(points),3)
    for i in 1:N
        cpoints[i,:] = Vector(- cog + points[i])
    end
    return (cog,cpoints)
end

function fitplane(points::Points)
    (cog,cpoints) = center(points)
    (U,S,V) = svd(cpoints)

    n = V[:,3]
    r = n*vecdot(cog,n)[1]

    @assert abs(norm(n)-1) < 1e-8
    @assert abs(abs(normalized(r-cog)⋅n)) < 1e-14
    @assert abs(abs(n⋅normalized(r))-1) < 1e-14
    @assert norm(r) < norm(cog)

    return Plane(n,r,points.A)
end

function fitline(points::Points)
    (cog,cpoints) = center(points)
    (U,S,V) = svd(cpoints)
    v = V[:,1]
    r = cog-vecdot(v,cog).*v

    @assert abs(norm(v)-1) < 1e-8
    @assert abs(abs(normalized(r-cog)⋅v) - 1) < 1e-14
    @assert abs(v⋅normalized(r)) < 1e-14
    @assert norm(r) < norm(cog)

    return Line(Point(v,points.A), Point(r,points.A),points.A)
end


function project(plane::Plane, points::Points)
    checkframes(plane,points)
    N = size(points)
    n = plane.n
    P = n.p*n.p'
    ppoints = zeros(points)
    for i = 1:N
        ppoints[i] = points[i] + (- P*points[i] + plane.r)
    end
    return ppoints
end

function project(plane::Plane, line::Line)
    checkframes(plane,line)
    n = plane.n
    n_projplane = normalized(n × line.v)
    P = n.p*n.p'
    vp = normalized(line.v - P*line.v)
    @assert abs(vp⋅n_projplane) < 1e-14

    A = [vp',n', n_projplane']
    b = [0, plane.r⋅n, line.r⋅n_projplane]
    rp = Point(A\b,line.A)
    @assert abs(rp ⋅ vp) < 1e-14
    @assert abs(vp ⋅ n) < 1e-14
    return Line(vp, rp)
end

function project(line::Line, points::Points)
    checkframes(points,line)
    N = size(points)
    v = line.v
    ppoints = zeros(points)
    for i = 1:N
        ppoints[i] = ((points[i]⋅v).*v  + line.r)
    end
    return ppoints
end



function framefromfeatures(feature1, feature2, origin, B)
    checkframes(feature1[2],feature2[2])
    x = y = z = false

    r1 = findfirst(lowercase(feature1[1][1]) .== ['x','y','z'])
    r2 = findfirst(lowercase(feature2[1][1]) .== ['x','y','z'])
    if r1 == 0 || r2 == 0; error("Your feature directions ($(feature1[1]), $(feature2[1])) are not recognized"); end
    r3 = 6-r1-r2
    s1 = feature1[1][2] == '+' ? 1 : -1
    s2 = feature2[1][2] == '+' ? 1 : -1

    R1 = Frames.normalized(Frames.getdir(feature1[2]).p * s1) |> Vector
    R2 = Frames.normalized(Frames.getdir(feature2[2]).p * s2) |> Vector
    R3 = Frames.normalized(R1×R2) |> Vector
    R2 = R3×R1

    _R = zeros(3,3)
    _R[:,r1] = R1
    _R[:,r2] = R2
    _R[:,r3] = R3

    (det(_R) < 0) &&  (_R[:,r3] *= -1)
    @assert abs(det(_R)-1) < 1e-14
    R = Mat(_R)

    Frame(R,origin.p,feature1[2].A, B)


end

function getdir(f::GeometricObject)
    if isa(f,Plane)
        return f.n
    elseif isa(f,Line)
        return f.v
    elseif isa(f,Frame)
        return Point(Rx(f),f.A)
    end
end



function readcloud(file)
    coordlines = prepfile(file)
    points = Points("T")
    for line in coordlines
        p = Point(line[2:4],"T")
        push!(points,p)
    end
    println("Read $(length(points)) points into cloud with reference T")
    points
end

function readplane(file)
    coordlines = prepfile(file)
    T = readTmatrix(coordlines)
    n = T[1:3,3]
    o = T[1:3,4]
    r = (n⋅o)*n
    Plane(n,r,"T")

end

function readTmatrix(coordlines)
    T = eye(4,4)
    for (i,line) in enumerate(coordlines)
        T[i,:] = line
    end
    T
end

function prepfile(file)
    f = open(file)
    text = readlines(f)
    close(f)
    coordlines = Array[]
    for line in text
        stripline = strip(line,[' ','[',']','"',';','\r','\n',' '])
        coords = split(stripline,[',',' '])
        coords = filter(x-> x!="",coords)
        try
            fcoords = float(coords)
            isempty(fcoords) && continue
            push!(coordlines,fcoords)
        catch e
            break
        end
    end
    coordlines
end


end
