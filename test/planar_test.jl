using LinearAlgebra
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    # Initial Conditions
    r0 = [0.64; 0.0; 0.76]
    v0 = [-0.48; 0.0; 0.0]
    th0 = deg2rad(60.0) # Remove q0 initial condition in subproblem to recreate plot
    q0 = [cos(th0 / 2); 0.0; sin(th0 / 2); 0.0]
    w0 = [0.0; 0.0; 0.0]
    m0 = 1
    x0 = [r0; v0; q0; w0; m0]

    # Rocket Mass properties
    g = 0.108
    r = 2 / 781.02
    h = 20 / 781.02
    diag = [(m0 / 12) * (3 * r^2 + h^2); (m0 / 12) * (3 * r^2 + h^2); 0.5 * m0 * r^2]
    Inertia = Diagonal(diag)
    Isp = 1 / (0.0738 * 9.807)
    cg = 10 / 781.02
    mdry = 0.277

    # Constraints
    Fmin = 0.024
    Fmax = 0.164
    gs = deg2rad(45.0)
    thmax = deg2rad(60.0)
    wmax = 143.84
    dmax = deg2rad(10.0)

    # Hyperparameters
    K = 50
    Nsub = 10
    wD = 1
    wDσ = 1
    wnu = 9e5

    # Set up problem
    par = PTR.PARAMS(x0, g, mdry, cg, Inertia, Isp, Fmin, Fmax, gs, thmax, wmax, dmax)
    p = PTR.ptr(K, par)
    p.Nsub = Nsub
    p.wD = wD
    p.wDσ = wDσ
    p.wnu = wnu

    # Solve problem
    PTR.solveTraj!(p)
    PTR.animatePlanarTrajectory(p)
    PTR.plotall(p)
end