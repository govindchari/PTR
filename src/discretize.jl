function FOH_discretize(p::ptr)

    for k = 1:p.K-1
        idx = p.idx
        println(k)

        # Get integrator initial conditions
        P0 = zeros(idx.N)
        P0[idx.x] = p.xref[:, k]
        P0[idx.phi] = reshape(I(p.nx), (p.nx^2, 1))
        P0[idx.Bm] = zeros(p.nx * p.nu)
        P0[idx.Bp] = zeros(p.nx * p.nu)
        P0[idx.S] = zeros(p.nx)
        P0[idx.z] = zeros(p.nx)

        # Integrate State
        z = RK4(df, P0, (k - 1) * p.dτ, p.dτ, p.Nsub, p)

        # Package discrete time matrices (Multiply by Ak)
        p.xprop[:, k] = z[idx.x]
        Ak = reshape(z[idx.phi], (p.nx, p.nx))
        p.A[:, :, k] = Ak
        p.Bm[:, :, k] = Ak * reshape(z[idx.Bm], (p.nx, p.nu))
        p.Bp[:, :, k] = Ak * reshape(z[idx.Bp], (p.nx, p.nu))
        p.S[:, k] = Ak * z[idx.S]
        p.z[:, k] = Ak * z[idx.z]

    end

end
function getState(τ::Float64, u::Function, p::ptr)
    # Uses RK4 to integrate propagate state from previous node to point between nodes
    k = Int(floor(τ / p.dτ)) + 1
    t0 = (k - 1) * p.dτ
    dt = τ - t0
    if (dt == 0)
        return p.xref[:, k]
    else
        df(t, x, p) = p.σref * p.f(x, u(t, p))
        xprop = RK4(df, p.xref[:, k], t0, dt, p.Nsub, p)
        return xprop
    end
end
function uprop(τ::Float64, p::ptr)
    k = Int(floor(τ / p.dτ)) + 1
    lm = (k * p.dτ - τ) / p.dτ
    lp = (τ - (k - 1) * p.dτ) / p.dτ
    if (k == p.K)
        return p.uref[:, k]
    else
        return lm * p.uref[:, k] + lp * p.uref[:, k+1]
    end
end
function df(τ::Float64, P::Array{Float64,1}, p::ptr)

    idx = p.idx
    k = Int(floor(τ / p.dτ)) + 1
    # Linearized Dynamics
    lm(t) = (k * p.dτ - t) / p.dτ
    lp(t) = (t - (k - 1) * p.dτ) / p.dτ
    # uprop(t) = lm(t) * p.uref[:, k] + lp(t) * p.uref[:, k+1]

    hfidx(t) = Int(t / (p.dτ / p.Nsub)) + 1
    A(t) = p.σref * p.dfx(getState(t, uprop, p), uprop(t, p))
    B(t) = p.σref * p.dfu(getState(t, uprop, p), uprop(t, p))
    S(t) = p.f(getState(t, uprop, p), uprop(t, p))
    z(t) = -reshape(A(t), (p.nx, p.nx)) * getState(t, uprop, p) - reshape(B(t), (p.nx, p.nu)) * uprop(t, p)


    Px = P[idx.x]
    Pphi = P[idx.phi]
    PBm = P[idx.Bm]
    PBp = P[idx.Bp]
    PS = P[idx.S]
    Pz = P[idx.z]

    Pphi_mat = reshape(Pphi, (p.nx, p.nx))
    F = factorize(Pphi_mat)

    dP = zeros(idx.N)
    dP[idx.x] = p.σref * p.f(Px, uprop(τ, p))
    dP[idx.phi] = reshape(A(τ) * Pphi_mat, (p.nx^2, 1))
    dP[idx.Bm] = reshape(F \ (lm(τ) * B(τ)), (p.nx * p.nu, 1))
    dP[idx.Bp] = reshape(F \ (lp(τ) * B(τ)), (p.nx * p.nu, 1))
    dP[idx.S] = F \ S(τ)
    dP[idx.z] = F \ z(τ)

    return dP

end
function RK4(dz::Function, z0::Array{Float64,1}, t0::Float64, dt::Float64, Nsub::Int64, p::ptr)
    h = dt / Nsub
    z = copy(z0)
    for i = 1:Nsub
        # hifi_idx = Int(t0 / dt) + 1 + i
        τ = t0 + (i - 1) * h
        k1 = dz(τ, z, p)
        k2 = dz(τ + 0.5 * h, z + 0.5 * h * k1, p)
        k3 = dz(τ + 0.5 * h, z + 0.5 * h * k2, p)
        k4 = dz(τ + h, z + h * k3, p)
        z += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        # p.xref_hifi[:, hifi_idx] = z[p.idx.x]
    end
    return z
end

