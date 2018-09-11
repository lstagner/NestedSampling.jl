using MultimodalNestedSampling
using SpecialFunctions

invphi(p) = sqrt(2) * erfinv(2 * p - 1.0)

function generate_data(α,β,n=100)
    θ = -pi/2 .+ pi*rand(n)
    x = α .+ β*tan.(θ)
    return x
end

n = 5
d = zeros(n)

d .= generate_data(0.5,2.0,n)

lighthouse_log_like_all = function loglike_all(cube)
    alpha = 2*invphi(cube[1])
    beta = 5*cube[2]

    l = sum(log(beta/pi) .- log.(beta^2 .+ (d .- alpha).^2))
    return l
end

problem = MultiNest(lighthouse_log_like_all, 2, "lighthouse")

run(problem,nlive=100)
