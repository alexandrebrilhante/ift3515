f(x) = x^2 + x - 2 * sqrt(x);

df(x) = 2 * x + 1 - 1 / sqrt(x);

function bisection(f::Function, a::Float64, b::Float64, δ::Float64 = 1e-6)
    k::Int64 = 1
    if a > b
        c = a
        a = b
        b = c
    end
    fa = f(a)
    fb = f(b)
    if fa * fb > 0
        println("The function must be of opposite signs at the bounds.")
        return
    end
    d = b - a
    c = a + d / 2
    fc = f(c)
    while d > δ
        k += 1
        if fc == 0
            a = b = c
            break
        elseif fc * fa < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
        d = b - a
        c = a + d / 2
        fc = f(c)
    end
    return k, fc, a, b
end

@show bisection(df, 0.0, 1.0)
