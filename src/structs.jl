struct IDX
    N::Int64
    x::UnitRange{Int64}
    phi::UnitRange{Int64}
    Bm::UnitRange{Int64}
    Bp::UnitRange{Int64}
    S::UnitRange{Int64}
    z::UnitRange{Int64}
    function IDX(nx, nu)
        idx_x = 1:nx
        idx_phi = (nx+1):(nx+nx^2)
        idx_Bm = (nx+nx^2+1):(nx+nx^2+nx*nu)
        idx_Bp = (nx+nx^2+nx*nu+1):(nx+nx^2+nx*nu+nx*nu)
        idx_S = (nx+nx^2+2*nx*nu+1):(nx+nx^2+2*nx*nu+nx)
        idx_z = (2*nx+nx^2+2*nx*nu+1):(2*nx+nx^2+2*nx*nu+nx)

        new(nx^2 + 3 * nx + 2 * nx * nu, idx_x, idx_phi, idx_Bm, idx_Bp, idx_S, idx_z)
    end

end
struct ptr

    # Problem paramters
    nx::Int64
    nu::Int64
    K::Int64
    dτ::Float64

    # Hyperparameters
    Nsub::Int64

    # Dynamics and Jacobians
    f::Function
    dfx::Function
    dfu::Function

    # Reference Trajectories
    xref::Array{Float64,2}
    uref::Array{Float64,2}
    σref::Float64

    # Discrete Dynamics
    idx::IDX
    xprop::Array{Float64,2}
    # xref_hifi::Array{Float64,2}
    A::Array{Float64,3}
    Bm::Array{Float64,3}
    Bp::Array{Float64,3}
    S::Array{Float64,2}
    z::Array{Float64,2}

    function ptr(nx, nu, K, Nsub, f, dfx, dfu)
        new(nx, nu, K, 1 / (K - 1), Nsub, f, dfx, dfu, zeros(nx, K), zeros(nu, K), 1.0, IDX(nx, nu), zeros(nx, K - 1), zeros(nx, nx, K - 1), zeros(nx, nu, K - 1), zeros(nx, nu, K - 1), zeros(nx, K - 1), zeros(nx, K - 1))
    end
end
