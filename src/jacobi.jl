using Optim

function jacobi(f::Function, x0::Vector, h::Float64, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    fsearch(α) = f(x-α*e[j]);
    while k < 500
        for j = 1:n
            e = zeros(n)
            e[j] = 1.0
            α = Optim.minimizer(optimize(fsearch, 0, h, GoldenSection()))
            x += α
        end
        k += 1
    end
    return x, k
end

@show jacobi(f, zeros(2), 2.0)
