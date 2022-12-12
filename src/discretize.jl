function FOH_discretize(p::ptr)
    A(τ) -> p.σref[τ] * p.dfx(p.xref[τ], p.uref[τ])
    B(τ) -> p.σref[τ] * p.dfu(p.xref[τ], p.uref[τ])
    S(τ) -> p.f(p.xref[τ], p.uref[τ])
    z(τ) -> -A(τ) * p.xref[τ] - B(τ) * p.uref[τ]

end
function RK4()

end
