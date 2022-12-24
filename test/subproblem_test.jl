using LinearAlgebra
using ForwardDiff
using Plots
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    g = 9.807
    g0 = 9.807
    Isp = 311.0
    gvec = [0; 0; -g]
    r_arm = [0; 0; -1]
    a = 1 / (Isp * g0)
    Inertia = I(3)

    # Rocket Dynamics (x = [r v q w m])
    function f(x, u)
        v = x[4:6]
        q = x[7:10]
        w = x[11:13]
        m = x[14]
        B = [-q[2] -q[3] -q[4]
            q[1] -q[4] -q[3]
            q[4] q[1] -q[2]
            -q[3] q[2] q[1]]
        dr = v
        dv = u / m + gvec
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

    nx = 14
    nu = 3
    K = 11
    Nsub = 10
    r0 = [0.0;0.0;1000.]
    v0 = [0.0;0.0;0.]
    q0 = [1.0;0.0;0.0;0.0]
    w0 = [0.0;0.0;0.0]
    m0 = 100
    x0 = [r0;v0;q0;w0;m0]

    p = PTR.ptr(nx, nu, K, f, dfx, dfu, x0)
    PTR.FOH_discretize!(p)
    PTR.initialize!(p)
    PTR.solveSubproblem!(p)

end