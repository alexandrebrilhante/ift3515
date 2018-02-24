using Optim

f(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2;

function gaussseidel(f::Function, x0::Vector, h::Float64, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    fsearch(α) = f(y - α * e[j]);
    while k < 500
        y = x
        for j = 1:n
            e = zeros(n)
            e[j] = 1.0
            α = Optim.minimizer(optimize(fsearch, 0, h, GoldenSection()))
            y += α * e
        end
        x = y
        k += 1
    end
    return x, k
end

@show gauss_seidel(f, zeros(2), 2.0)
