using LinearAlgebra
using ForwardDiff
using Plots
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    # g = 9.807
    # g0 = 9.807
    # Isp = 311.0
    # gvec = [0; 0; -g]
    # r_arm = [0; 0; -1]
    # a = 1 / (Isp * g0)
    # Inertia = I(3)

    nx = 14
    nu = 3
    K = 15
    Nsub = 10
    r0 = [0.64; 0.0; 0.76]
    v0 = [-0.48; 0.0; 0.0]
    q0 = [sqrt(2)/2; 0.0; 0.0; sqrt(2)/2]
    w0 = [0.0; 0.0; 0.0]
    m0 = 1
    x0 = [r0; v0; q0; w0; m0]

    g = 0.108
    gvec = [0; 0; -g]
    r_arm = [0; 0; -1e-3]
    a = 0.0738
    r = 6 / 781.02
    h = 30 / 781.02
    diag = [0.5 * m0 * r^2; (m0 / 12) * (3 * r^2 + h^2); (m0 / 12) * (3 * r^2 + h^2)]
    Inertia = Diagonal(diag)

    # Rocket Dynamics (x = [r v q w m])
    function f(x, u)
        v = x[4:6]
        q = x[7:10]
        q0 = q[1]
        q1 = q[2]
        q2 = q[3]
        q3 = q[4]
        w = x[11:13]
        m = x[14]
        B = [-q[2] -q[3] -q[4]
            q[1] -q[4] -q[3]
            q[4] q[1] -q[2]
            -q[3] q[2] q[1]]
        bCi = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2)
            2*(q1*q2-q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3+q0*q1)
            2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) q0^2-q1^2-q2^2+q3^2]
        dr = v
        dv = bCi' * u / m + gvec
        dq = 0.5 * B * w
        dw = Inertia \ (cross(r_arm, u) - cross(w, I * w))
        dm = -a * norm(u)
        return [dr; dv; dq; dw; dm]
    end
    function dfx(x, u)
        return ForwardDiff.jacobian(dx -> f(dx, u), x)
    end
    function dfu(x, u)
        return ForwardDiff.jacobian(du -> f(x, du), u)
    end
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

    p = PTR.ptr(nx, nu, K, f, dfx, dfu, x0)
    PTR.solveTraj!(p)
    plot(p.xref[1, :], p.xref[3, :])
    quiver!(p.xref[1, :], p.xref[3, :], quiver=(p.uref[1, :], p.uref[3, :]))

    un = []
    tvc = []
    for k=1:K
        append!(un, norm(p.uref[:,k]))
        append!(tvc, acos(p.uref[3,k]/norm(p.uref[:,k])))
    end
    plot(rad2deg.(tvc))
end