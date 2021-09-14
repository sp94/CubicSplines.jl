using Test
using CubicSplines
using Random

@testset "cubic_splines.jl" begin

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

    # Interpolate a sine function with white noise
    xdata = range(0,stop=4π,length=20) .+ 0.5*rand(MersenneTwister(1234),20)
    ydata = sin.(xdata)

    # Create the spline
    spline = CubicSpline(xdata, ydata)
    @test spline[π] ≈ 0.00652422 rtol = 1e-5
    @test spline[0.35] ≈ 0.346190 rtol = 1e-5

    # Non-allocating vector filling
    x = [1.0, 2.0]
    y = zeros(2)
    spline(y,x)
    for i in 1:length(y)
        @test y[i] == spline(x[i])
    end

    # Plot the spline
    # xs = range(xdata[1], stop=xdata[end], length=1000)
    # ys = spline[xs]
    # plt.plot(xdata, ydata, "o")
    # plt.plot(xs, ys)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.close()

    # Extrapolate the spline with a linear curve to the left
    spline = CubicSpline(xdata, ydata; extrapl=[1.0]) # Linear extrapolation with an incline of 1.0

    # Data from outside the interpolation range
    left_value = -1.0
    right_value = 14.0
    @test spline[left_value] ≈ -1.004278 rtol = 1e-5 # Extrapolate linearly to the left
    @test_throws String spline[right_value] # No extrapolation (throw an error)

    # Plot the spline
    # xs = range(left_value, stop=xdata[end], length=1000)
    # ys = spline[xs]
    # plt.plot(xdata, ydata, "o")
    # plt.plot(xs, ys)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.title("With linear extrapolation to the left")
    # plt.close()

    # Extrapolate the spline with a quadratic curve to the left and a cubic curve to the right
    spline = CubicSpline(xdata, ydata; extrapl=[0.0; 0.5], extrapr=[1.0; 0.0; 0.5])
    left_value = -1.0
    right_value = 14.0
    @test spline[left_value] ≈ 1.13020 rtol = 1e-5 # Extrapolate linearly to the left
    @test spline[right_value] ≈ 1.93991 rtol = 1e-5 # Extrapolate linearly to the left

    # Plot the spline
    # xs = range(left_value, stop=right_value, length=2000)
    # ys = spline[xs]
    # plt.plot(xdata, ydata, "o")
    # plt.plot(xs, ys)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.title("With extrapolation on both sides")
    # plt.close()

    # Gradient calculation
    xdata = range(0,stop=4pi,length=20)
    ydata = sin.(xdata)
    spline = CubicSpline(xdata, ydata; extrapl=[2.0; 1.0], extrapr=[1.0; 0.0; 0.5])

    # Plot the spline
    # xs = range(left_value, stop=right_value, length=2000)
    # ys = spline[xs]
    # plt.plot(xdata, ydata, "o")
    # plt.plot(xs, ys)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.title("With extrapolation on both sides")
    # plt.close()

    # Calculate the spline gradient at different points and evaluate it
    @test gradient(spline, -2.0) ≈ -2.0
    @test gradient(spline, -2.0, 2) ≈ 2.0
    @test gradient(spline, -2.0, 3) ≈ 0.0
    @test gradient(spline, pi/2, 1) ≈ 0.01851987 rtol = 1e-5
end
