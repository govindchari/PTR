function plotall(p::ptr)
    un = []
    tvc = []
    th = []
    qnorm = []
    u_norm = []
    ui = zeros(p.nu, p.K)
    for k = 1:p.K
        append!(un, norm(p.uref[:, k]))
        append!(tvc, acos(p.uref[3, k] / norm(p.uref[:, k])))
        q0 = p.xref[p.idx.q, k][1]
        q1 = p.xref[p.idx.q, k][2]
        q2 = p.xref[p.idx.q, k][3]
        q3 = p.xref[p.idx.q, k][4]
        bCi = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q0*q3) 2*(q1*q3-q0*q2)
            2*(q1*q2-q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3+q0*q1)
            2*(q1*q3+q0*q2) 2*(q2*q3-q0*q1) q0^2-q1^2-q2^2+q3^2]
        ui[:, k] = bCi' * p.uref[:, k]
        append!(th, 2 * acos(clamp(q0, -1, 1)))
        append!(qnorm, norm(p.xref[p.idx.q, k]))
        append!(u_norm, norm(p.uref[:, k]))
    end
    pygui(true)

    plt.figure(1)
    ax = plt.axes(projection="3d")
    ax.plot3D(p.xref[1, :], p.xref[2, :], p.xref[3, :], color="black", label="Trajectory")
    ax.quiver(p.xref[1, :], p.xref[2, :], p.xref[3, :], ui[1, :], ui[2, :], ui[3, :], length=0.5, color="r", label="Thrust Vector")
    ax.set_xlabel("Downrange")
    ax.set_ylabel("Crossrange")
    ax.set_zlabel("Altitude")
    ax.set_title("PDG Trajectory")

    plt.figure(2)
    plt.plot(rad2deg.(th), color="black", label="Alignment Angle")
    plt.axhline(y=rad2deg(p.par.thmax), color="r", linestyle="--", label="Constraint")
    plt.title("Alignment")
    plt.xlabel("Time")
    plt.ylabel("Angle (deg)")
    plt.grid()

    plt.figure(3)
    plt.plot(u_norm, color="black", label="Thrust Magnitude")
    plt.axhline(y=p.par.Fmin, color="r", linestyle="--", label="Constraint")
    plt.axhline(y=p.par.Fmax, color="r", linestyle="--", label="Constraint")
    plt.title("Thrust Magnitude")
    plt.xlabel("Time")
    plt.ylabel("Thrust")
    plt.grid()

    plt.figure(4)
    plt.plot(rad2deg.(tvc), color="black", label="Gimbal Angle")
    plt.axhline(y=rad2deg(p.par.dmax), color="r", linestyle="--", label="Constraint")
    plt.title("Gimbal Angle")
    plt.xlabel("Time")
    plt.ylabel("Angle (deg)")
    plt.grid()

    plt.legend()
    plt.show()

end