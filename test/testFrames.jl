using Robotlib.Frames
R = toOrthoNormal(randn(3,3))
t = randn(3)
f = Frame{1}(R,t)
F = f
