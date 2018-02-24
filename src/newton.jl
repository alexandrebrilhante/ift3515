using ForwardDiff

f(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2;

g(x) = ForwardDiff.gradient(f, x);

H(x) = ForwardDiff.hessian(f, x);

function newton(f::Function, g::Function, H::Function,
                x0::Vector, δ::Float64 = 1e-6)
    k::Int64 = 0
    x = x0
    n::Int6 = length(x)
    dfx = 1.0
    while norm(dfx) > δ && k < 500
        dfx = g(x)
        d2fx = H(x)
        x -= d2fx \ dfx
        k += 1
    end
    return x, k
end

@show newton(f, g, H, zeros(2))
