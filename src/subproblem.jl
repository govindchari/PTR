function solveSubproblem!(p::ptr)
    x = Variable(p.nx, p.K)
    u = Variable(p.nu, p.K)
    σ = Variable(1)

    nu = Variable(p.nx * (p.K - 1))
    Δ = Variable(p.K)
    Δσ = Variable(1)

    dx = x - p.xref
    du = u - p.uref
    dσ = σ - p.σref

    idx = p.idx
    par = p.par
    r = x[idx.r, :]
    v = x[idx.v, :]
    q = x[idx.q, :]
    w = x[idx.w, :]
    m = x[idx.m, :]

    # Objective
    objective = (σ + p.wD * norm(Δ) + p.wDσ * norm(Δσ, 1) + p.wnu * sumsquares(nu)) / p.wnu

    # Dynamics Constraint
    constraints = Constraint[
        x[:, k+1] == p.A[:, :, k] * x[:, k] + p.Bm[:, :, k] * u[:, k] + p.Bp[:, :, k] * u[:, k+1] + p.S[:, k] * σ + p.z[:, k] + nu[(k-1)*p.nx+1:k*p.nx] for k in 1:p.K-1
    ]

    # Boundary Conditions
    # x = [r v q w m]
    push!(constraints, x[idx.r, 1] == par.x0[idx.r])
    push!(constraints, x[idx.v, 1] == par.x0[idx.v])
    push!(constraints, x[idx.q, 1] == par.x0[idx.q])
    push!(constraints, x[idx.w, 1] == par.x0[idx.w])
    push!(constraints, x[idx.m, 1] == par.x0[idx.m])

    push!(constraints, x[idx.r, p.K] == zeros(3))
    push!(constraints, x[idx.v, p.K] == zeros(3))
    push!(constraints, x[idx.q, p.K] == [1.0; 0.0; 0.0; 0.0])
    push!(constraints, x[idx.w, p.K] == zeros(3))

    # State Constraints
    for k = 1:p.K
        push!(constraints, par.mdry - m[k] <= 0)
        push!(constraints, norm([1 0 0; 0 1 0] * r[:, k]) * tan(par.gs) - [0 0 1] * r[:, k] <= 0)
        push!(constraints, norm(w[:, k]) <= par.wmax)
        push!(constraints, cos(par.thmax) <= 1 - 2 * (sumsquares([0 1 0 0; 0 0 1 0] * q[:, k])))
    end

    # Control Constraints
    Xi(k) = p.uref[:, k] / norm(p.uref[:, k])
    for k = 1:p.K
        push!(constraints, norm(u[:, k]) <= par.Fmax)
        push!(constraints, par.Fmin <= dot(Xi(k), u[:, k]))
        push!(constraints, cos(par.dmax) * norm(u[:, k]) <= u[3, k])
    end

    # Trust Regions
    for k = 1:p.K
        push!(constraints, norm(dot(dx[:, k], dx[:, k]) + dot(du[:, k], du[:, k]), 2) <= Δ[k])
    end
    push!(constraints, norm(dσ, 2) <= Δσ)

    prob = minimize(objective, constraints)
    solve!(prob, ECOS.Optimizer, silent_solver=true)
    p.xref = evaluate(x)
    p.uref = evaluate(u)
    p.σref = evaluate(σ)
    p.vc = reshape(evaluate(nu), (p.nx, p.K - 1))
    p.Δ = evaluate(Δ)
    p.Δσ = evaluate(Δσ)
end