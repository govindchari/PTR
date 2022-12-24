function initialize!(p::ptr)
    idx = p.idx
    g = 9.807

    p.xref[:, 1] = p.x0
    p.xref[:, p.K] = p.xT
    p.xref[idx.m, p.K] = p.mdry
    q0 = p.x0[idx.q]
    qf = p.xT[idx.q]
    th = acos(dot(q0, qf))

    # Initialize State and Control
    for k = 2:p.K-1
        p.xref[idx.r, k] = p.x0[idx.r] + ((p.xT[idx.r] - p.x0[idx.r]) / (p.K - 1)) * (k - 1)
        p.xref[idx.v, k] = p.x0[idx.v] + ((p.xT[idx.v] - p.x0[idx.v]) / (p.K - 1)) * (k - 1)
        if (abs(th) > 1e-6)
            p.xref[idx.q, k] = (q0 * sin((1 - (k - 1) * p.dτ) * th) + qf * sin((k - 1) * p.dτ * th)) / sin(th)
        else
            p.xref[idx.q, k] = q0
        end
        p.xref[idx.w, k] = p.x0[idx.w] + ((p.xT[idx.w] - p.x0[idx.w]) / (p.K - 1)) * (k - 1)
        p.xref[idx.m, k] = p.x0[idx.m] + ((p.mdry - p.x0[idx.m]) / (p.K - 1)) * (k - 1)
        p.uref[:, k] = [0.0; 0.0; p.xref[idx.m, k] * g] # (warning) gravity is hardcoded
    end
    p.σref += sqrt(2 * norm(p.x0[idx.r]) / g)
end

function solveTraj(p::ptr)
end