## Example usage

```
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