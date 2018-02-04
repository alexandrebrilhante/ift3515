f(x) = x^2 + x - 2 * sqrt(x);

N::Int64 = 50
F = ones(N)

for i = 3:N
    F[i] = F[i-1] + F[i-2]
end

function fibonacci(f::Function, a::Float64, b::Float64)
    k::Int64 = 1
    d = b - a
    xL = a + (F[N-2] / F[N]) * d
    xR = a + (F[N-1] / F[N]) * d
    fL = f(xL)
    fR = f(xR)
    while k < N - 2
        k += 1
        if fL < fR
            b = xR
            d = b - a
            xR = xL
            fR = fL
            xL = a + (F[N-k-1] / F[N-k+1]) * d
            fL = f(xL)
        elseif fL > fR
            b = xL
            d = b - a
            xL = xR
            fL = fR
            xR = a + (F[N-k] / F[N-k+1]) * d
            fR = f(xR)
        elseif fL == fR
            k += 1
            if k < N - 2
                a = xL
                b = xR
                d = b - a
                xL = a + (F[N-k-1] / F[N-k+1]) * d
                xR = a + (F[N-k] / F[N-k+1]) * d
                fL = f(xL)
                fR = f(xR)
            end
        end
    end
    return a, b
end

@show fibonacci(f, 0.0, 1.0)
