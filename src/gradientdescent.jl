using ForwardDiff, Optim

f(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(f, x);

function gradientdescent(f::Function, g::Function, x0::Vector,
                         h::Float64, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    dfx = ones(n)
    fsearch(α) = f(x - α * dfx);
    while norm(dfx) > δ && k < 500
        α = Optim.minimizer(optimize(fsearch, 0, h, GoldenSection()))
        dfx = g(x)
        x -= α * dfx
        k += 1
    end
    return x, k
end

@show gradientdescent(f, g, zeros(2), 2.0)
