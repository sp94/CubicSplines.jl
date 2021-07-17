#=
Samuel J Palmer
Implementation is based on:
[1]: Akima, Hiroshi (1970). "A new method of interpolation and smooth curve fitting
    based on local procedures". Journal of the ACM. 17: 589–602.

Original module by Samuel J. Palmer, see https://github.com/sp94/CubicSplines.jl
Modified by Stefan Mathis to allow for customized extrapolation and to use a
binary search for interval identification.
=#

module CubicSplines

export CubicSpline, gradient

# As described in section 2.3 of 1
function extrapolate(xs, ys, extrap)

    # Get first / last three values
    x1, x2, x3 = xs
    y1, y2, y3 = ys

    # Extrapolate the x-coordinates of two additional datapoints with (8) from [1]
    x4 = x3-x1+x2
    x5 = x3-x1+x3

    # Quadratic extrapolation of (x4,y4), (x5,y5)
    if extrap === nothing
        y4 = y3 + 2*(x4-x3)*(y3-y2)/(x3-x2) - (x4-x3)*(y2-y1)/(x2-x1)
        y5 = y4 + 2*(x5-x4)*(y4-y3)/(x4-x3) - (x5-x4)*(y3-y2)/(x3-x2)
    else
        # Use the evalpoly macro to calculate the extrapolation values
        y4 = @evalpoly(x4-x3, y3, extrap...)
        y5 = @evalpoly(x5-x3, y3, extrap...)
    end
    return x4, x5, y4, y5
end

"""
Construct an Akima spline according to [1]. In addition to the original method,
this implementation also allows for customized extrapolation at both spline
sides. The extrapolation behaviour is defined by the coefficients of a 3rd degree
polynomial which can be given as additional keyword arguments extrapl and extrapr
to the constructor. extrapl and/or extrapr should be either "nothing" or a
set of coefficients.

Case 1: Extrapolation is set to "nothing"
================================================================================
In this case, the default quadratic extrapolation according to section 2.3 of [1]
is used to construct the spline and the spline throws an error when evaluated
outside its boundaries.

Case 2: Extrapolation is set to a vector
================================================================================
The container is given as a n-element collection of values:
    extrapl = [p1, p2, p3, ..., pn]
which are used to evaluate the following function:
    y = p0 + p1*(x-x3) + p2*(x-x3)^2 + p3*(x-x3)^3 + ... + pn*(x-x3)^n
The constant p0 is equal to y3 by definition (continuous spline). x3 and y3
are the first/last datapoint values respectively.
"""
struct CubicSpline{T <: Real}
    xs::Vector{T}
    ys::Vector{T}
    ts::Vector{T}
    ps::Matrix{T} # Coefficients of the polynom for each spline section
    logx::Bool
    extrapl::Union{Nothing,Vector{T}}
    extrapr::Union{Nothing,Vector{T}}
    function CubicSpline(xs, ys; logx=false, extrapl=nothing, extrapr=nothing)

        # Check input data
        @assert length(xs) == length(ys)
        @assert issorted(xs)
        if !issorted(xs)
            idx = sortperm(xs)
            xs = xs[idx]
            ys = ys[idx]
        end
        if logx
            xs = log.(xs)
        end
        if length(xs) < 5
            throw("At least five datapoints are necessary to construct a cubic Akima spline, however only $(length(xs)) datapoints were given.")
        end

        xb, xa, yb, ya = extrapolate(xs[3:-1:1], ys[3:-1:1], extrapl)
        xy, xz, yy, yz = extrapolate(xs[end-2:end], ys[end-2:end], extrapr)
        xs = [xa; xb; xs; xy; xz]
        ys = [ya; yb; ys; yy; yz]
        ms = zeros(eltype(ys), length(xs)-1)
        for i in 1:length(ms)
            ms[i] = (ys[i+1]-ys[i]) ./ (xs[i+1]-xs[i])
        end
        ts = zeros(eltype(ys), length(xs))
        for i in 3:length(xs)-2

            # Equals (1) in [1]
            m1, m2, m3, m4 = ms[i-2:i+1]
            # As described in [1] p.591 (parentheses block), this is an arbitrary
            # convention to guarantee uniqueness of the solution
            if m1 == m2 && m3 == m4
                ts[i] = (m2+m3)/2
            else
                numer = abs(m4-m3)*m2 + abs(m2-m1)*m3
                denom = abs(m4-m3)    + abs(m2-m1)
                ts[i] = numer ./ denom
            end
        end

        # Spline coefficients
        ps = zeros(eltype(ys), length(xs)-1, 4)
        for i in 3:length(xs)-2
            x1, x2 = xs[i:i+1]
            y1, y2 = ys[i:i+1]
            t1, t2 = ts[i:i+1]
            p0 = y1
            p1 = t1
            p2 = (3(y2-y1)/(x2-x1)-2t1-t2)/(x2-x1)
            p3 = (t1+t2-2(y2-y1)/(x2-x1))/(x2-x1)^2
            ps[i,:] = [p0, p1, p2, p3]
        end

        # Create extrapolation coefficients vector
        if extrapl !== nothing
            extrapl = vcat([ys[3]],[c for c in extrapl]) # 3rd value equals first value of the original dataset
        end
        if extrapr !== nothing
            extrapr = vcat([ys[end-2]],[c for c in extrapr]) # 3rd last value equals last value of the original dataset
        end
        return new{eltype(xs)}(xs[3:end-2], ys[3:end-2], ts[3:end-2], ps[3:end-2,:], logx, extrapl, extrapr)
    end
end
# Functor to make the CubicSpline class callable
function (obj::CubicSpline)(x)
    return obj[x]
end
# Non-allocating version for vectors
function (obj::CubicSpline)(y::AbstractVector,x::AbstractVector)
    for i in 1:length(x)
        y[i] = obj[x[i]]
    end
    return y
end

# As described in section 2.2 of [1]
function Base.getindex(spline::CubicSpline, x::Real)

    if spline.logx
        x = log(x)
    end
    if isnan(x)
        return NaN
    end

    # Value sits on upper interval border
    if x == spline.xs[end]
        idx = length(spline.xs)-1

    # Find the corresponding interval
    else

        # Find the interval index
        idx = _binary_search_interval(spline.xs, x)

        # Extrapolation to the "left"
        if idx == 0
            if spline.extrapl === nothing
                throw("x too small ($x < $(spline.xs[1]))")
            else
                return @inbounds @evalpoly(x-spline.xs[1], spline.extrapl...)
            end
        end

        # Extrapolation to the "right"
        if idx == length(spline.xs)
            if spline.extrapr === nothing
                throw("x too big ($x > $(spline.xs[end]))")
            else
                return @inbounds @evalpoly(x-spline.xs[end], spline.extrapr...)
            end
        end
    end

    # Evaluate the spline segment
    return @inbounds @evalpoly(x-spline.xs[idx], spline.ps[idx,:]...)
end

function Base.getindex(spline::CubicSpline, xs::AbstractArray{<:Real,1})
    return [spline[x] for x in xs]
end

function Base.firstindex(spline::CubicSpline)
    return spline.xs[1]
end

function Base.lastindex(spline::CubicSpline)
    return spline.xs[end]
end

"""
Evaluate the spline gradient of degree one at a given point x
"""
function gradient(spline::CubicSpline, x::Real)
    return gradient(spline, x, 1)
end

"""
Generalized gradient of degree n at a given point x
"""
function gradient(spline::CubicSpline, x::Real, n::Integer)

    # Value sits on upper interval border
    if x == spline.xs[end]
        idx = length(spline.xs)-1

    # Find the corresponding interval
    else

        # Find the interval index
        idx = _binary_search_interval(spline.xs, x)

        # Extrapolation to the "left"
        if idx == 0
            if spline.extrapl === nothing
                throw("x too small ($x < $(spline.xs[1]))")
            else
                coeff = spline.extrapl
                x_ref = spline.xs[1]
            end

        # Extrapolation to the "right"
        elseif idx == length(spline.xs)
            if spline.extrapr === nothing
                throw("x too big ($x > $(spline.xs[end]))")
            else
                coeff = spline.extrapr
                x_ref = spline.xs[end]
            end
        else
            coeff = spline.ps[idx,:]
            x_ref = spline.xs[idx]

        end
    end

    # Calculate the coefficients of degree n
    coeff = _diffpoly(coeff, n, 1)

    # If the coefficients are empty, the derivative is 0.
    if isempty(coeff)
        return zero(x)
    else
        return @evalpoly(x-x_ref, coeff...)
    end
end

"""
Binary search to find the interval which contains a given number. The interval
is given as a sorted array of real numbers. Each neighbouring pair is interpreted
as an interval. The function returns the lower index of the interval containing
the given number.

If the returned index is 0, the number is located in the interval (-Inf, array[1])
If the returned index is length(array), the number is located in the interval [array[end], Inf)
"""
@inline function _binary_search_interval(array::Vector{T}, number::Real) where T <: Real

    left = 1
    right = length(array)

    @inbounds while left <= right
        center = (left+right)÷2
        if array[center] > number
            right = center - 1
        else # array[center] <= number
            left = center + 1
        end
    end
    return left - 1
end

"""
Calculate the coefficients for the differentiated polygon
"""
function _diffpoly(coeff::AbstractArray{<:Real,1}, n::Integer=1, counter::Integer=1)

    # Discard the first coefficient, since it vanishes during the differentiation process
    coeff = coeff[2:end] .* collect(1:(length(coeff)-1))

    # If the differentiation degree is larger than the counter, differentiate the
    # differentiated coefficients again.
    if n == counter
        return coeff
    else
        return _diffpoly(coeff, n, counter+1)
    end
end

end # module
