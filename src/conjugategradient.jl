A = [4 2; 2 4]

b = [0; 1]

x = [1; 1]

function conjugategradient(A::Matrix, b::Vector, x0::Vector, δ::Float64 = 1e-6)
    n::Int64 = length(x0)
    x = x0
    g = b + A * x
    d = -g
    k::Int64 = 0
    while norm(g) > δ
        Ad = A * d
        normd = dot(d, Ad)
        α = -dot(d, g) / normd
        x += α * d
        g = b + A * x
        β = dot(g, Ad) / normd
        d = -g + β * d
        k += 1
    end
    normd = dot(d, A * d)
    α = -dot(d, g) / normd
    x += α * d
    return x, k
end

@show conjugategradient(A, b, x)
