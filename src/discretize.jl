function FOH_discretize(p::ptr)

    for k = 1:p.K-1
        idx = p.idx

        # Get integrator initial conditions
        P0 = zeros(idx.N)
        P0[idx.x] = p.xref[:, k]
        P0[idx.phi] = reshape(I(nx), (nx^2, 1))
        P0[idx.Bm] = zeros(p.nx * p.nu)
        P0[idx.Bp] = zeros(p.nx * p.nu)
        P0[idx.S] = zeros(p.nx)
        P0[idx.z] = zeros(p.nx)

        # Integrate State
        z = RK4(df, P0, k*dτ, dτ, p.Nsub, p)

        # Package discrete time matrices (Multiply by Ak)
    end

end
function df(τ::Float64, P::Array{Float64,1}, p::ptr)

    # Linearized Dynamics
    A(t) = p.σref[t] * p.dfx(p.xref[t], p.uref[t])
    B(t) = p.σref[t] * p.dfu(p.xref[t], p.uref[t])
    S(t) = p.f(p.xref[t], p.uref[t])
    z(t) = -A(t) * p.xref[t] - B(t) * p.uref[t]
    lm(_k, t) = ((_k + 1) * p.dτ - t) / p.dτ
    lp(_k, t) = (t - _k * p.dτ) / p.dτ
    uprop(_k, t) = lm(_k, t) * p.uref[:, _k] + lp(_k, t) * p.uref[:, _k+1]

    Px = P[idx.x]
    Pphi = P[idx.phi]
    PBm = P[idx.Bm]
    PBp = P[idx.Bp]
    PS = P[idx.S]
    Pz = P[idx.z]

    Pphi_mat = reshape(Pphi, (nx, nx))
    F = LU(Pphi_mat)

    dP = zeros(idx.N)
    dP[idx.x] = p.f(Px, uprop(τ), p.σref)
    dP[idx.phi] = reshape(A(τ) * Pphi_mat, (p.nx^2, 1))
    dP[idx.Bm] = reshape(F \ (lm(k, τ) * B(τ)), (p.nx*p.nu, 1))
    dP[idx.Bp] = reshape(F \ (lp(k, τ) * B(τ)), (p.nx*p.nu, 1))
    dP[idx.S] = F \ S(τ)
    dP[idx.z] = F \ z(τ)

    return dP

end
function RK4(dz::Function, z0::Array{Float64,1}, t0::Float64, dt::Float64, Nsub::Int64, p::ptr)
    h = dt/Nsub
    z = copy(z0)
    for i=1:Nsub
        τ = t0 + (i-1)*h
        k1 = dz(τ,z,p)
        k2 = dz(τ + 0.5*h, z + 0.5*h*k1,p)
        k3 = dz(τ + 0.5*h, z + 0.5*h*k2,p)
        k4 = dz(τ + h, z + h*k3,p)
        z += (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    return z
end

