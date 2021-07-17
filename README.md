# Cubic Splines

A simple package for interpolating 1D data with Akima cubic splines, based on "A New Method of Interpolation and Smooth Curve Fitting Based on Local Parameters", Akima, 1970.

Works for both uniformly and non-uniformly spaced data points.

## Example usage

```julia
using PyPlot
using CubicSplines

xdata = range(0,stop=4pi,length=20) .+ 0.5rand(20)
ydata = sin.(xdata)
plot(xdata, ydata, "o")

spline = CubicSpline(xdata, ydata)

xs = range(xdata[1], stop=xdata[end], length=1000)
ys = spline[xs]
plot(xs, ys)
xlabel("x")
ylabel("y")
```

![Example sinusoid](img/example_sinusoid.png)

## Extrapolation

It is possible to add an individual extrapolating polynomial on both sides of the spline. When
constructing the spline, a smooth transition between spline and polynomial is achieved
by enforcing equal values as well as first and second derivatives at `xdata[1]` and
`xdata[end]` for the spline and the extrapolating polynomial. For each of the
two spline ends, an individual polynomial can be given (or none at all, resulting
in an error if values from that region are requested). The left-hand side
polynomial coefficients are given as keyword argument `extrapl`, the right-hand
side coefficients as `extrapr` (assuming a x-axis where -Inf is on the left
and +Inf on the right).

Examples for extrapolation from runtests.jl:

```julia
using PyPlot
using CubicSplines
using Test

# Extrapolate the spline with a quadratic curve to the left and a cubic curve to the right
spline = CubicSpline(xdata, ydata; extrapl=[0.0; 0.5], extrapr=[1.0; 0.0; 0.5])
left_value = -1.0
right_value = 14.0
@test spline[left_value] ≈ 1.13020 rtol = 1e-5 # Extrapolate linearly to the left
@test spline[right_value] ≈ 1.93991 rtol = 1e-5 # Extrapolate linearly to the left

# Plot the spline
plt.figure()
xs = range(left_value, stop=right_value, length=2000)
ys = spline[xs]
plt.plot(xdata, ydata, "o")
plt.plot(xs, ys)
plt.xlabel("x")
plt.ylabel("y")
plt.title("With extrapolation on both sides")

# Extrapolate the spline with a linear curve to the left
spline = CubicSpline(xdata, ydata; extrapl=[1.0], extrapr=nothing) # Linear extrapolation with an incline of 1.0

# Data from outside the interpolation range
left_value = -1.0
right_value = 14.0
@test spline[left_value] ≈ -1.004278 rtol = 1e-5 # Extrapolate linearly to the left
@test_throws String spline[right_value] # No extrapolation (throw an error)

# Plot the spline
plt.figure()
xs = range(left_value, stop=xdata[end], length=1000)
ys = spline[xs]
plt.plot(xdata, ydata, "o")
plt.plot(xs, ys)
plt.xlabel("x")
plt.ylabel("y")
plt.title("With linear extrapolation to the left")
```

The extrapolation polynomials are defined by their coefficients and the first / last
value of `xdata` respectively:

```julia
extrapl = [p1, p2, p3, ..., pn]
y = p0 + p1*(x-xdata[1]) + p2*(x-xdata[1])^2 + p3*(x-xdata[1])^3 + ... + pn*(x-xdata[1])^n
```

## Gradient calculation

The gradient of grade n at value x can be requested by `gradient(spline, x, n)`.
This also works for the extrapolation polynomials
