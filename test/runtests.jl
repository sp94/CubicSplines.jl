using Test

using CubicSplines

# Create data that is symmetric about zero
xdata = -2.0:1.0:2.0
ydata = mod.(xdata, 2)
@assert reverse(xdata) == -xdata
@assert reverse(ydata) == +ydata

# Test spline is symmetric about zero
xs = -2.0:0.5:2.0
ys = CubicSpline(xdata, ydata)[xs]
@assert reverse(xs) == -xs
@test reverse(ys) == +ys
