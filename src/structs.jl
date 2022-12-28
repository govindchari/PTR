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
    Fmin::Float64  # Minimum throttle
    Fmax::Float64  # Max throttle
    gs::Float64    # Glideslope
    thmax::Float64 # Max deviation from vertical
    wmax::Float64  # Max angular rate
    dmax::Float64  # Max gimbal angle

    function PARAMS(x0, g, mdry, Fmin, Fmax, gs, thmax, wmax, dmax)
        new(x0, g, mdry, Fmin, Fmax, gs, thmax, wmax, dmax)
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

    function ptr(nx::Int64, nu::Int64, K::Int64, f::Function, dfx::Function, dfu::Function, par::PARAMS)
        Nsub = 10
        wD = 1
        wDσ = 1
        wnu = 1e4
        new(nx, nu, K, 1 / (K - 1), Nsub, wD, wDσ, wnu, f, dfx, dfu, zeros(nx, K), zeros(nu, K), 0.0, zeros(nx, K), IDX(nx, nu), zeros(nx, K - 1), zeros(nx, nx, K - 1), zeros(nx, nu, K - 1), zeros(nx, nu, K - 1), zeros(nx, K - 1), zeros(nx, K - 1), par)
    end
end
