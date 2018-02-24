using ForwardDiff

f(x::Vector) = -10 * x[1]^2 + 10 * x[2]^2 + 4 * sin(x[1] * x[2]) - 2 * x[1] + x[1]^4;

g(x) = ForwardDiff.gradient(f, x);

l = [0.0; 0.0]

u = [10.0; 10.0]

function armijo(f::Function, dfx::Vector, x::Vector, d::Vector,
                α::Float64 = 1.0, β1::Float64 = 1e-4, β2::Float64 = 0.9)
    s = β1 * dot(dfx, d)
    fx = f(x)
    fxcand = f(x + α * d)
    while fxcand > fx + α * s
        α *= β2
        fxcand = f(x + α * d)
    end
    return α
end

function projection!(x::Vector, l::Vector, u::Vector)
    n = length(x)
    for i = 1:n
        x[i] = max(min(u[i], x[i]), l[i])
    end
    return x
end

function projectedgradient(f::Function, g::Function, x0::Vector,
                           l::Vector, u::Vector, δ::Float64 = 1e-6,
                           maxiter::Int64 = 1000)
    k::Int64 = 0
    x = x0
    n::Int64 = length(x)
    d = ones(n)
    dfx = ones(n)
    while norm(d) > δ && k < maxiter
        dfx = g(x)
        d = projection!(x - dfx, l, u) - x
        α = armijo(f, dfx, x, d)
        x += α * d
        k += 1
    end
    @show x, k
end

@show projectedgradient(f, g, zeros(2), l, u)
