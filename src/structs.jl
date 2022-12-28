struct IDX
    N::Int64
    x::UnitRange{Int64}
    phi::UnitRange{Int64}
    Bm::UnitRange{Int64}
    Bp::UnitRange{Int64}
    S::UnitRange{Int64}
    z::UnitRange{Int64}
    r::UnitRange{Int64}
    v::UnitRange{Int64}
    q::UnitRange{Int64}
    w::UnitRange{Int64}
    m::Int64

    function IDX(nx, nu)
        idx_x = 1:nx
        idx_phi = (nx+1):(nx+nx^2)
        idx_Bm = (nx+nx^2+1):(nx+nx^2+nx*nu)
        idx_Bp = (nx+nx^2+nx*nu+1):(nx+nx^2+nx*nu+nx*nu)
        idx_S = (nx+nx^2+2*nx*nu+1):(nx+nx^2+2*nx*nu+nx)
        idx_z = (2*nx+nx^2+2*nx*nu+1):(2*nx+nx^2+2*nx*nu+nx)

        new(nx^2 + 3 * nx + 2 * nx * nu, idx_x, idx_phi, idx_Bm, idx_Bp, idx_S, idx_z, 1:3, 4:6, 7:10, 11:13, 14)
    end

end
struct PARAMS
    # Initial Condition
    x0::Array{Float64,1}
    g::Float64 # Gravitatinal constraints

    # Rocket parameters/constraints
    mdry::Float64  # Dry mass
    cg::Float64
    Inertia::Diagonal{Float64,Vector{Float64}}
    Isp::Float64
    Fmin::Float64  # Minimum throttle
    Fmax::Float64  # Max throttle
    gs::Float64    # Glideslope
    thmax::Float64 # Max deviation from vertical
    wmax::Float64  # Max angular rate
    dmax::Float64  # Max gimbal angle

    function PARAMS(x0, g, mdry, cg, Inertia, Isp, Fmin, Fmax, gs, thmax, wmax, dmax)
        new(x0, g, mdry, cg, Inertia, Isp, Fmin, Fmax, gs, thmax, wmax, dmax)
    end
end
mutable struct ptr

    # Problem paramters
    nx::Int64
    nu::Int64
    K::Int64
    dτ::Float64

    # PTR Hyperparameters
    Nsub::Int64
    wD::Float64
    wDσ::Float64
    wnu::Float64

    # Dynamics and Jacobians
    f::Function
    dfx::Function
    dfu::Function

    # Reference Trajectories
    xref::Array{Float64,2}
    uref::Array{Float64,2}
    σref::Float64
    vc::Array{Float64,2}
    Δ::Array{Float64, 1}
    Δσ::Float64

    # Discrete Dynamics
    idx::IDX
    xprop::Array{Float64,2}
    A::Array{Float64,3}
    Bm::Array{Float64,3}
    Bp::Array{Float64,3}
    S::Array{Float64,2}
    z::Array{Float64,2}

    # Problem paramters
    par::PARAMS

    function ptr(K::Int64, par::PARAMS)
        nx = 14
        nu = 3
        Nsub = 10
        wD = 1
        wDσ = 1
        wnu = 5e4

        g0 = 9.807
        a = 1 / (par.Isp * g0)
        r_arm = [0; 0; -par.cg]

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
            dv = bCi' * u / m + [0; 0; -par.g]
            dq = 0.5 * B * w
            dw = par.Inertia \ (cross(r_arm, u) - cross(w, par.Inertia * w))
            dm = -a * norm(u)
            return [dr; dv; dq; dw; dm]
        end
        function dfx(x, u)
            return ForwardDiff.jacobian(dx -> f(dx, u), x)
        end
        function dfu(x, u)
            return ForwardDiff.jacobian(du -> f(x, du), u)
        end

        new(nx, nu, K, 1 / (K - 1), Nsub, wD, wDσ, wnu, f, dfx, dfu, zeros(nx, K), zeros(nu, K), 0.0, zeros(nx, K), zeros(K), 0.0, IDX(nx, nu), zeros(nx, K - 1), zeros(nx, nx, K - 1), zeros(nx, nu, K - 1), zeros(nx, nu, K - 1), zeros(nx, K - 1), zeros(nx, K - 1), par)
    end
end
