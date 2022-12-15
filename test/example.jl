using LinearAlgebra
using Plots

include("../src/PTR.jl")
using .PTR

let
    function f(x::Array{Float64,1}, u::Array{Float64,1})
        return [x[2]; -sin(x[1]) + u[1]]
    end
    function fint(t::Float64, x::Array{Float64,1}, p::ptr)
        return [x[2]; -sin(x[1])]
    end
    function dfx(x::Array{Float64,1}, u)
        return [0 1;-cos(x[1]) 0]
    end
    function dfu(x::Array{Float64,1}, u)
        return [0;1]
    end
    nx = 2
    nu = 1
    K = 11
    Nsub = 10
    p = ptr(nx,nu,K,Nsub,f,dfx,dfu)
    println(p.dτ)
    FOH_discretize(p)

    # dt = 10.0
    # Nsub = 100
    # z0 = [1.0;0.0]
    # T = 100
    # t0 = 0.0
    # z_list = zeros(nx,Int(T/dt))
    # t_list = 1.0:Int(T/dt)
    # z_list[:,1] = z0
    # for i = 2:Int(T/dt)
    #     z_list[:,i] = RK4(fint, z_list[:,i-1], t_list[i], dt, Nsub, p)
    # end
    # print(z_list)
    # plot(z_list[1,:])

end