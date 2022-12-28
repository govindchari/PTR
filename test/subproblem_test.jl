using LinearAlgebra
using ForwardDiff
using Plots
using BenchmarkTools

include("../src/PTR.jl")
using .PTR

let
    K = 15
    r0 = [0.64; 0.0; 0.76]
    v0 = [-0.48; 0.0; 0.0]
    th0 = deg2rad(55.0)
    q0 = [cos(th0 / 2); 0.0; 0.0; sin(th0 / 2)]
    w0 = [0.0; 0.0; 0.0]
    m0 = 1
    x0 = [r0; v0; q0; w0; m0]

    g = 0.108
    r = 2 / 781.02
    h = 20 / 781.02
    diag = [0.5 * m0 * r^2; (m0 / 12) * (3 * r^2 + h^2); (m0 / 12) * (3 * r^2 + h^2)]
    Inertia = Diagonal(diag)
    Isp = 1 / (0.0738 * 9.807)
    cg = 1e-3
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
    un = []
    tvc = []
    th = []
    qnorm = []
    nu_norm = []
    ui = zeros(p.nu, p.K)
    for k = 1:K
        append!(un, norm(p.uref[:, k]))
        append!(tvc, acos(p.uref[3, k] / norm(p.uref[:, k])))
        q0 = p.xref[p.idx.q, k][1]
        q1 = p.xref[p.idx.q, k][2]
        q2 = p.xref[p.idx.q, k][3]
        q3 = p.xref[p.idx.q, k][4]
        bCi = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2)
            2*(q1*q2-q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3+q0*q1)
            2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) q0^2-q1^2-q2^2+q3^2]
        ui[:, k] = bCi' * p.uref[:, k]
        append!(th, 2 * acos(clamp(q0, -1, 1)))
        append!(qnorm, norm(p.xref[p.idx.q, k]))
        if k != K
            append!(nu_norm, norm(p.vc[:, k]))
        end
    end
    plot(p.xref[1, :], p.xref[3, :])
    quiver!(p.xref[1, :], p.xref[3, :], quiver=(ui[1, :], ui[3, :]))
    # plot(qnorm)
    plot(nu_norm)
    # plot(rad2deg.(th))
    # plot!(ones(p.K)*rad2deg(p.par.thmax))
end