# Samuel J Palmer
# Based on "A New Method of Interpolation and Smooth Curve Fitting Based on Local Parameters", Akima, 1970

module CubicSplines

function extrapolate(xs, ys)
    # Quadratic extrapolation of (x4,y4), (x5,y5)
    x1, x2, x3 = xs
    y1, y2, y3 = ys
    x4 = x3-x1+x2
    x5 = x3-x1+x3
    y4 = y3 + 2*(x4-x3)*(y3-y2)/(x3-x2) - (x4-x3)*(y2-y1)/(x2-x1)
    y5 = y4 + 2*(x5-x4)*(y4-y3)/(x4-x3) - (x5-x4)*(y3-y2)/(x3-x2)
    return x4, x5, y4, y5
end

struct CubicSpline
    xs
    ys
    ts
    logx
    function CubicSpline(xs, ys; logx=false)
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
        xb, xa, yb, ya = extrapolate(xs[3:-1:1], ys[3:-1:1])
        xy, xz, yy, yz = extrapolate(xs[end-2:end], ys[end-2:end])
        xs = [xa; xb; xs; xy; xz]
        ys = [ya; yb; ys; yy; yz]
        ms = zeros(eltype(ys), length(xs)-1)
        for i in 1:length(ms)-1
            ms[i] = (ys[i+1]-ys[i]) ./ (xs[i+1]-xs[i])
        end
        ts = zeros(eltype(ys), length(xs))
        for i in 3:length(xs)-2
            m1, m2, m3, m4 = ms[i-2:i+1]
            if m1 == m2 && m3 == m4
                ts[i] = (m2+m3)/2
            else
                numer = abs(m4-m3)*m2 + abs(m2-m1)*m3
                denom = abs(m4-m3)    + abs(m2-m1)
                ts[i] = numer ./ denom
            end
        end
        return new(xs[3:end-2], ys[3:end-2], ts[3:end-2], logx)
    end
end

function Base.getindex(spline::CubicSpline, x::Real)
    if spline.logx
        x = log(x)
    end
    if x < spline.xs[1]
        throw("x too small ($x < $(spline.xs[1]))")
    elseif x > spline.xs[end]
        throw("x too big ($x > $(spline.xs[end]))")
    end
    i = 1
    while !(spline.xs[i] <= x <= spline.xs[i+1])
        i += 1
    end
    x1, x2 = spline.xs[i:i+1]
    y1, y2 = spline.ys[i:i+1]
    t1, t2 = spline.ts[i:i+1]
    p0 = y1
    p1 = t1
    p2 = (3(y2-y1)/(x2-x1)-2t1-t2)/(x2-x1)
    p3 = (t1+t2-2(y2-y1)/(x2-x1))/(x2-x1)^2
    return p0 + p1*(x-x1) + p2*(x-x1)^2 + p3*(x-x1)^3
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

export CubicSpline

end # module
