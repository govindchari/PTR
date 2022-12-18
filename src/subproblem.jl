function solveSubproblem(p::ptr)
    x = Variable(p.nx, p.K)
    u = Variable(p.nu, p.K)
    σ = Variable(1)

    v = Variable(p.nx, p.K - 1)
    D = Variable(p.K)
    Dσ = Variable(1)
    
    dx = x - p.xref
    du = u - p.uref
    dσ = σ - p.σref

    # Boundary Conditions

    # Dynamics Constraint
    constraints = Constraint[
        x[:, k+1] == p.A[:, :, k] * x[:, k] + p.Bm[:, :, k] * u[:, k] + p.Bp[:, :, k] * u[:, k+1] + p.S[:, k] * σ + z[:, k] + v[:, k] for k in 1:p.K-1
    ]

    # State Constraints

    # Control Constraints

    # Trust Regions


end