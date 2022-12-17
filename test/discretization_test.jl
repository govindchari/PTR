using LinearAlgebra
using Plots
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    # Pendulum Dynamics and Jacobians
    w = 2 * pi
    function f(x::Array{Float64,1}, u::Array{Float64,1})
        return [x[2]; -w^2 * sin(x[1]) + u[1]]
    end
    function dfx(x::Array{Float64,1}, u)
        return [0 1; -w^2*cos(x[1]) 0]
    end
    function dfu(x::Array{Float64,1}, u)
        return [0; 1]
    end

    nx = 2
    nu = 1
    K = 11
    Nsub = 10
    p = PTR.ptr(nx, nu, K, Nsub, f, dfx, dfu)
    @btime PTR.FOH_discretize($p)

    x0 = [0.1; 0.0]
    u = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]'
    function df(τ::Float64, z::Array{Float64,1}, p::PTR.ptr)
        # Function for integrator
        k = Int(floor(τ / p.dτ)) + 1
        lm = (k * p.dτ - τ) / p.dτ
        lp = (τ - (k - 1) * p.dτ) / p.dτ
        if (k == p.K)
            u_ = u[:, k]
        else
            u_ = lm * u[:, k] + lp * u[:, k+1]
        end
        return f(z, u_)
    end

    xd = zeros(nx, K)
    xc = zeros(nx, K)
    xd[:, 1] = x0
    xc[:, 1] = x0
    z = x0[1] * cos.((w) * (0:0.1:1))
    for i = 1:K-1
        xd[:, i+1] = p.A[:, :, i] * xd[:, i] + p.Bm[:, :, i] * u[:, i] + p.Bp[:, :, i] * u[:, i+1]
        xc[:, i+1] = PTR.RK4(df, xc[:, i], (i - 1) * p.dτ, p.dτ, 1, p)
    end

    plot(0:10, xd[1, :])
    plot!(0:10, xc[1, :])
    # println(norm(xd-xc, Inf))
    # plot!(0:10, z)

end