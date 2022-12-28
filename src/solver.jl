function initialize!(p::ptr)
    idx = p.idx
    par = p.par

    p.xref[:, 1] = par.x0
    p.xref[idx.r, p.K] = zeros(3)
    p.xref[idx.v, p.K] = zeros(3)
    p.xref[idx.q, p.K] = [1.0; 0.0; 0.0; 0.0]
    p.xref[idx.w, p.K] = zeros(3)

    p.xref[idx.m, p.K] = par.mdry

    p.uref[:, 1] = [0.0; 0.0; p.xref[idx.m, 1] * par.g]
    p.uref[:, p.K] = [0.0; 0.0; p.xref[idx.m, p.K] * par.g]
    q0 = par.x0[idx.q]
    qf = p.xref[idx.q, p.K]
    th = acos(dot(q0, qf))


    # Initialize State and Control
    for k = 2:p.K-1
        p.xref[idx.r, k] = par.x0[idx.r] + ((p.xref[idx.r, p.K] - par.x0[idx.r]) / (p.K - 1)) * (k - 1)
        p.xref[idx.v, k] = par.x0[idx.v] + ((p.xref[idx.v, p.K] - par.x0[idx.v]) / (p.K - 1)) * (k - 1)

        # Slerp
        if (abs(th) > 1e-6)
            p.xref[idx.q, k] = (q0 * sin((1 - (k - 1) * p.dτ) * th) + qf * sin((k - 1) * p.dτ * th)) / sin(th)
        else
            p.xref[idx.q, k] = q0
        end

        p.xref[idx.w, k] = par.x0[idx.w] + ((p.xref[idx.w, p.K] - par.x0[idx.w]) / (p.K - 1)) * (k - 1)
        p.xref[idx.m, k] = par.x0[idx.m] + ((par.mdry - par.x0[idx.m]) / (p.K - 1)) * (k - 1)
        p.uref[:, k] = [0.0; 0.0; p.xref[idx.m, k] * par.g]
    end
    p.σref += sqrt(2 * norm(par.x0[idx.r]) / par.g)
end

function solveTraj!(p::ptr)
    initialize!(p)
    for _ = 1:10
        FOH_discretize!(p)
        solveSubproblem!(p)
    end
end