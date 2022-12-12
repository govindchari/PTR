mutable struct ptr
    nx::Int64
    nu::Int64
    K::Int64
    dτ::Float64

    # Dynamics and Jacobians
    f
    dfx
    dfu

    xref::Array{Float64,2}
    uref::Array{Float64,2}
    σref::Array{Float64,1}

    xprop::Array{Float64,2}
    A::Array{Float64,3}
    Bm::Array{Float64,3}
    Bp::Array{Float64,3}
    S::Array{Float64,2}
    z::Array{Float64,2}

    function ptr(nx, nu, K, f, dfx, dfu)
        new(nx, nu, K, 1 / (K - 1), f, dfx, dfu, zeros(nx, K), zeros(nu, K), zeros(K), zeros(K), zeros(nx, nx, K), zeros(nx, nu, K), zeros(nx, nu, K), zeros(nx, K), zeros(nx, K))
    end
end
