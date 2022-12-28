using LinearAlgebra
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    K = 15
    r0 = [0.64; 0.0; 0.76]
    v0 = [-0.48; 0.0; 0.0]
    th0 = deg2rad(60.0)
    q0 = [cos(th0 / 2); 0.0; sin(th0 / 2); 0.0]
    w0 = [0.0; 0.0; 0.0]
    m0 = 1
    x0 = [r0; v0; q0; w0; m0]

    g = 0.108
    r = 2 / 781.02
    h = 20 / 781.02
    diag = [0.5 * m0 * r^2; (m0 / 12) * (3 * r^2 + h^2); (m0 / 12) * (3 * r^2 + h^2)]
    Inertia = Diagonal(diag)
    Isp = 1 / (0.0738 * 9.807)
    cg = 10 / 781.02
    mdry = 0.277
    Fmin = 0.024
    Fmax = 0.164
    gs = deg2rad(45.0)
    thmax = deg2rad(90.0)
    wmax = 143.84
    dmax = deg2rad(10.0)

    par = PTR.PARAMS(x0, g, mdry, cg, Inertia, Isp, Fmin, Fmax, gs, thmax, wmax, dmax)
    p = PTR.ptr(K, par)

    PTR.solveTraj!(p)
    PTR.plotall(p)
end