"""
Usage:
```julia
path = Pkg.dir("Robotlib","src","applications","frames")


add_frame_name!(\"SEAM\",\"Weld seam frame\")
add_frame_name!(\"TAB\",\"Table frame\")

T_RB_Tm = MAT.matread(path*\"T_RB_T.mat\")[\"T_RB_T\"]
T_TF_TCPm = MAT.matread(path*\"T_TF_TCP.mat\")[\"T_TF_TCP\"]
T_T_TABm = MAT.matread(path*\"T_T_Table.mat\")[\"T_T_Table\"]

T_RB_T = Frame(T_RB_Tm,\"RB\",\"T\")
T_S_D = Frame(T_TF_TCPm,\"S\",\"D\")
T_T_TAB = Frame(T_T_TABm,\"T\",\"TAB\")

cloud_seam = readcloud(path*\"CloudSeam_edge.txt\")
plane_seam = readplane(path*\"PlaneSeam_edge.txt\")
cloud_seam_projected = project(plane_seam,cloud_seam)
line_seam = fitline(cloud_seam_projected)

T_T_SEAM = framefromfeatures((\"z+\",line_seam),(\"y-\",plane_seam),cloud_seam_projected[1],\"SEAM\")
T_RB_SEAM = T_RB_T*T_T_SEAM
T_RB_TAB = T_RB_T*T_T_TAB
T_TAB_SEAM = inv(T_T_TAB)*T_T_SEAM

MAT.matwrite(path*\"T_TAB_SEAM.mat\",[\"T_TAB_SEAM\" => T_TAB_SEAM.T])
MAT.matwrite(path*\"T_T_SEAM.mat\",[\"T_T_SEAM\" => T_T_SEAM.T])
MAT.matwrite(path*\"T_RB_TAB.mat\",[\"T_RB_TAB\" => T_RB_TAB.T])
println(\"Wrote T_TAB_SEAM, T_T_SEAM, T_RB_TAB to files in \$path\")

cloud_seam_RB = T_RB_T*cloud_seam
cloud_seam_projected_RB = T_RB_T*cloud_seam_projected
plane_seam_RB = T_RB_T*plane_seam
line_seam_RB = T_RB_T*line_seam



plotframe(Frame(eye(4),\"RB\",\"U\"),200, label=true)

plotpoints!(cloud_seam_RB)
plotpoints!(cloud_seam_projected_RB)
plotline!(line_seam_RB,500,label=\"Line seam\")
plotplane!(plane_seam_RB,200,label=\"Plane seam\")
plotframe!(T_RB_SEAM,200, label=true)
plotframe!(T_RB_TAB,200, label=true)

xlabel!(\"x\")
ylabel!(\"y\")
zlabel!(\"z\")
```
"""
module Frames

using Plots
export Frame, Point, Plane, Points, Line, GeometricObject, add_frame_name!
export readcloud, readTmatrix, readplane, fitline, fitplane, framefromfeatures, project
export plotline, plotplane, plotpoints, plot3Dsmart, plotframe
export inv, *,+,-,/,\,transpose,ctranspose

import Base: det, print, zeros, length, size, getindex, setindex!, convert, push!, show, start, next, done, +, *, ⋅, .*, /, ./, -, ×, transpose, ctranspose, \, inv
# using LaTeXStrings




# Type definitions --------------------------------------------
# -------------------------------------------------------------
abstract GeometricObject

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

function convert(::Type{Array{Float64,2}}, p::Points)
    N = length(p)
    points = zeros(N,3)
    for i = 1:N
        points[i,:] = p[i]'
    end
    return points
end

function zeros(points::Points)
    p = Point()
    N = length(points)
    newpoints = [p for i in 1:N]
    Points(newpoints,points.A)
end

type Frame <: GeometricObject

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

end

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

checkframes(x::GeometricObject, y::GeometricObject) = x.A != y.A &&  error("The two geometric objects do not have the same reference frame ($x -> $x.A, $y -> $y.A)")

function show(f::Frame) showcompact(round(f.T[1:3,1:4],4)) ; show(string(frame_names[f.A])*"->"*string(frame_names[f.B])) end
print(f::Frame) = show(f)

type Plane <: GeometricObject
    n::Point
    r::Point
    A::String
    Plane(n,r,A="U") = new(n,r,A)
    Plane(n::Vector,r::Vector,A="U") = new(Point(n,A),Point(r,A),A)
    Plane(points::Points) = fitplane(points)
end

type Line <: GeometricObject
    v::Point
    r::Point
    A::String
    Line(v,r,A="U") = new(v,r,A)
    Line(points::Points) = fitline(points)
end




# Plot functions ----------------------------------------------
# -------------------------------------------------------------

plot3Dsmart!(x::Array{Float64,2};kw...) = plot!(x[:,1],x[:,2],x[:,3];kw...)

function plotframe!(f::Frame = Frame(), length=1.0; label=false)
    #Plots a frame using XYZ-RGB
    o = F2t(f)
    x = Rx(f)
    y = Ry(f)
    z = Rz(f)
    plot3Dsmart!([o o+x*length]',c=:red)
    plot3Dsmart!([o o+y*length]',c=:green)
    plot3Dsmart!([o o+z*length]',c=:blue)
    if label && f.A != ""
        po = o-length/4*(x+y+z)
        annotate!(po[1],po[2],po[3],print(f))
    end
end


plotframe!(f::Matrix, length=1.0; label=false) = plotframe!(Frame(f),length,label=label)

function plotline!(l::Line, length=1.0; label="")
    coords = [(l.r-length*l.v).p, (l.r+length*l.v).p]
    scatter!(l.r[1],l.r[2],l.r[3])
    plot3Dsmart!(coords)
    if label != ""
        po = l.r.p-length/4*normalized(l.r.p)
        annotate!(po[1],po[2],po[3],label)
    end
end


function plotplane!(p::Plane, length=1.0; label="")
    coords = [(p.r).p (p.r+length*p.n).p]'
    scatter!(p.r[1],p.r[2],p.r[3],m=:^)
    plot3Dsmart!(coords,linespec)
    if label != ""
        po = p.r.p-length/4*p.n.p
        annotate!(po[1],po[2],po[3],label)
    end
end


plotpoints!(p::Points;kw...) = plot3Dsmart!(convert(Array{Float64,2},p);kw...)

for fun in ["plot3Dsmart","plotframe","plotline","plotplane"]
    @eval begin
        function $(Symbol(fun))(args...;kwargs...)
            fig = plot()
            $(Symbol(fun,!))(args...;kwargs...)
            return fig
        end
    end
end

# Frame functions  --------------------------------------------
# -------------------------------------------------------------

*(p::Point, c::Real) = Point(c.*p.p, p.A)
*(c::Real, p::Point) = *(p, c)
.*(p::Point, c::Real) = *(p, c)
.*(c::Real, p::Point) = *(p, c)
*(p::Point, R::Matrix) = error("Post multiplication of ::Point with ::Matrix not supported")
*(R::Matrix, p::Point) = Point(R*p.p, p.A)
/(p::Point, c::Real) = Point(p.p./c, p.A)
/(c::Real,p::Point) = /(p,c)
./(p::Point, c::Real) = Point(p,c)
./(c::Real,p::Point) = /(p,c)
+(p1::Point, p2::Point) = Point(p1.p + p2.p, p1.A)
-(p1::Point, p2::Point) = Point(p1.p - p2.p, p1.A)
+(p1::Vector{Float64}, p2::Point) = p1 + p2.p
-(p1::Vector{Float64}, p2::Point) = p1 - p2.p
+(p2::Point, p1::Vector{Float64}) = Point(p1 + p2.p, p2.A)
-(p2::Point, p1::Vector{Float64}) = Point(p1 - p2.p, p2.A)
⋅(p1::Point, p2::Point) = p1.p ⋅ p2.p
×(p1::Point, p2::Point) = Point(p1.p × p2.p, p1.A)

transpose(p::Point) = p.p'
ctranspose(p::Point) = p.p'
convert(::Type{Vector{Float64}}, p::Point) = p.p

*(x::Frame, y::Frame) = Frame(x.T*y.T,x.A,y.B)
function *(f::Frame, p::Point)
    f.B != p.A &&  error("Reference frames does not match between $f and $p ($(f.B) != $(p.A))")
    pv = (f.T*[p.p;1])[1:3]
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
normalized(p::Point) = Point(p.p/norm(p.p), p.A)
function normalize!(p::Point) p.p/=norm(p.p); end

function ⊥(R::Matrix)
    (U,S,V) = svd(R)
    a = diagm([1, 1, sign(det(U*V'))])
    return U*a*V'
end

function ⊥(F::Frame)
    F2  = deepcopy(F)
    R   = ⊥(F2R(F2))
    F2.T[1:3,1:3] = R
    return F2
end


# Projections and fitting  ------------------------------------
# -------------------------------------------------------------

function center(points::Points)
    N = size(points)
    cog = zeros(3)
    for i = 1:N
        cog += points[i]
    end
    cog /= N

    cpoints = zeros(length(points),3)
    for i in 1:N
        cpoints[i,:] = (- cog + points[i])'
    end
    return (cog,cpoints)
end

function fitplane(points::Points)
    (cog,cpoints) = center(points)
    (U,S,V) = svd(cpoints)

    n = V[:,3]
    r = n*(cog'*n)[1]

    assert(abs(norm(n)-1) < 1e-8)
    assert(abs(abs(normalized(r-cog)⋅n)) < 1e-14)
    assert(abs(abs(n⋅normalized(r))-1) < 1e-14)
    assert(norm(r) < norm(cog))

    return Plane(n,r,points.A)
end

function fitline(points::Points)
    (cog,cpoints) = center(points)
    (U,S,V) = svd(cpoints)
    v = V[:,1]
    r = cog-(v'*cog).*v

    assert(abs(norm(v)-1) < 1e-8)
    assert(abs(abs(normalized(r-cog)⋅v) - 1) < 1e-14)
    assert(abs(v⋅normalized(r)) < 1e-14)
    assert(norm(r) < norm(cog))

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
    assert(abs(vp⋅n_projplane) < 1e-14)

    A = [vp',n', n_projplane']
    b = [0, plane.r⋅n, line.r⋅n_projplane]
    rp = Point(A\b,line.A)
    assert(abs(rp ⋅ vp) < 1e-14)
    assert(abs(vp ⋅ n) < 1e-14)
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
    R = zeros(3,3)

    r1 = findfirst(lowercase(feature1[1][1]) .== ['x','y','z'])
    r2 = findfirst(lowercase(feature2[1][1]) .== ['x','y','z'])
    if r1 == 0 || r2 == 0; error("Your feature directions ($(feature1[1]), $(feature2[1])) are not recognized"); end
    r3 = 6-r1-r2
    s1 = feature1[1][2] == '+' ? 1 : -1
    s2 = feature2[1][2] == '+' ? 1 : -1

    R1 = normalized(getdir(feature1[2]).p * s1)
    R2 = normalized(getdir(feature2[2]).p * s2)
    R3 = normalized(R1×R2)
    R2 = R3×R1

    R[:,r1] = R1
    R[:,r2] = R2
    R[:,r3] = R3
    (det(R) < 0) &&  (R[:,r3] *= -1)
    assert(abs(det(R)-1) < 1e-14)

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
