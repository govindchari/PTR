function solveSubproblem(p::ptr)
    x = Variable(p.nx, p.K)
    u = Variable(p.nu, p.K)
    σ = Variable(1)

    nu = Variable(p.nx * (p.K - 1))
    D = Variable(p.K)
    Dσ = Variable(1)

    dx = x - p.xref
    du = u - p.uref
    dσ = σ - p.σref

    r = x[1:3, :]
    v = x[4:6, :]
    q = x[7:10, :]
    w = x[11:13, :]
    m = x[14, :]

    # Objective
    objective = σ + p.wD * norm(D) + p.wDσ * norm(Dσ, 1) + p.wnu * norm(nu, 1)

    # Dynamics Constraint
    constraints = Constraint[
        x[:, k+1] == p.A[:, :, k] * x[:, k] + p.Bm[:, :, k] * u[:, k] + p.Bp[:, :, k] * u[:, k+1] + p.S[:, k] * σ + p.z[:, k] + nu[(k-1)*p.nx+1:k*p.nx] for k in 1:p.K-1
    ]

    # Boundary Conditions
    # x = [r v q w m]
    push!(constraints, x[:, 1] == p.x0)
    push!(constraints, x[:, p.K] == zeros(p.nx))
    push!(constraints, x[7, p.K] == 1.0)

    # State Constraints
    for k = 1:p.K
        push!(constraints, p.mdry - m[k] <= 0)
        push!(constraints, norm([0 1 0; 0 0 1] * r[:, k]) * tan(p.gs) - [1 0 0] * r[:, k] <= 0)
        # push!(constraints, cos(p.thmax) <= 1 - 2 * (q[2, k]^2 + q[3, k]^2))
        push!(constraints, cos(p.thmax) <= 1 - 2 * (sumsquares([0 1 0 0;0 0 1 0]*q[:,k])))
        push!(constraints, norm(w[:, k]) <= p.wmax)
    end

    # Control Constraints
    Xi(k) = p.uref[:, k] / norm(p.uref[:, k])
    for k = 1:p.K
        q1 = q[1, k]
        q2 = q[2, k]
        q3 = q[3, k]
        q4 = q[4, k]
        bhat = [sumsquares([1 0 0 0;0 1 0 0]*q[:,k]) - sumsquares([0 0 1 0;0 0 0 1]*q[:,k]); 2 * (q2 * q3 + q1 * q4); 2 * (q2 * q4 - q1 * q2)]
        push!(constraints, norm(u[:, k]) <= p.Fmax)
        push!(constraints, p.Fmin <= dot(Xi(k), u[:, k]))
        # push!(constraints, cos(p.dmax) * norm(u[:, k]) <= dot(bhat, u[:, k]))
    end

    # Trust Regions
    for k = 1:p.K
        push!(constraints, norm(dot(dx[:, k], dx[:, k]) + dot(du[:, k], du[:, k]), 2) <= D[k])
    end
    push!(constraints, norm(dσ, 2) <= Dσ)

    prob = minimize(objective, constraints)
    solve!(prob, ECOS.Optimizer)

end