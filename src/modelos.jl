# Arquivo contendo modelos para simulação
using StaticArrays
using DifferentialEquations
using RigidBodyDynamics

urdf_kuka = "kuka.urdf"
mecanismo_kuka = parse_urdf(Float64, urdf_kuka)
state_kuka = MechanismState(mecanismo_kuka)


"""
Robô serial com 2 graus de liberdade na vertical
"""
function robot2dof(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::Vector{Float64}) where {T<:AbstractMatrix, Y<: AbstractFloat}

    function myrobot(du, u, p, t)
        m = SVector{2}([23.902, 1.285])
        l= SVector{2}([.45, .45])
        r= SVector{2}([0.091, 0.048])
        I= SVector{2}([1.266, 0.093])
        g = 9.81
        θ = SVector{2}(u[1:2])
        dθ = SVector{2}(u[3:4])
        cos1 = cos(θ[1])
        cos2 = cos(θ[2])
        sin1 = sin(θ[1])
        sin2 = sin(θ[2])
        cos12 = cos(θ[1]+θ[2])

        #Matriz de inércia
        m11 = m[1] * r[1]^2 + m[2]* (l[1]^2 + r[2]^2 + 2*l[1] * r[2] * cos2)+ I[1] + I[2]
        m12 = m[2] * (r[2]^2 + l[1] *r[2] * cos2) + I[2]
        m22 = m[2] * r[2]^2 + I[2]
        M = SMatrix{2,2}([m11 m12; m12 m22])

        #Matriz C
        h = -m[2] * l[1] * r[2] * sin2
        c11 = h * dθ[2]
        c12 = h * (dθ[1] + dθ[2])
        c21 = -h * dθ[1]
        c22 = 0
        C = SMatrix{2,2}([c11 c12; c21 c22])

        #Matriz gravidade
        g2 = m[2] * r[2] * g * cos12
        g1 = (m[1] * r[1] + m[2] * l[1]) * g * cos1 + g2
        G = SVector{2}([g1; g2])

        e = SVector{2}(xr - θ)
        tau = kp*e - kv*dθ
        du[1:2] = dθ
        du[3:4] = inv(M)*(tau - C * dθ - G)
        du[5:6] = tau
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(2,sol2)
end


"""
Robô serial com 2 graus de liberdade na vertical com dados
"""
function robot2dof(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::Vector{Z}, vr::Vector{Z}) where {T<:AbstractMatrix, Y<: AbstractFloat, Z<:Function}

    function myrobot(du, u, p, t)
        m = SVector{2}([23.902, 1.285])
        l= SVector{2}([.45, .45])
        r= SVector{2}([0.091, 0.048])
        I= SVector{2}([1.266, 0.093])
        g = 9.81
        θ = SVector{2}(u[1:2])
        dθ = SVector{2}(u[3:4])
        cos1 = cos(θ[1])
        cos2 = cos(θ[2])
        sin1 = sin(θ[1])
        sin2 = sin(θ[2])
        cos12 = cos(θ[1]+θ[2])

        #Matriz de inércia
        m11 = m[1] * r[1]^2 + m[2]* (l[1]^2 + r[2]^2 + 2*l[1] * r[2] * cos2)+ I[1] + I[2]
        m12 = m[2] * (r[2]^2 + l[1] *r[2] * cos2) + I[2]
        m22 = m[2] * r[2]^2 + I[2]
        M = SMatrix{2,2}([m11 m12; m12 m22])

        #Matriz C
        h = -m[2] * l[1] * r[2] * sin2
        c11 = h * dθ[2]
        c12 = h * (dθ[1] + dθ[2])
        c21 = -h * dθ[1]
        c22 = 0
        C = SMatrix{2,2}([c11 c12; c21 c22])

        #Matriz gravidade
        g2 = m[2] * r[2] * g * cos12
        g1 = (m[1] * r[1] + m[2] * l[1]) * g * cos1 + g2
        G = SVector{2}([g1; g2])

        e = SVector{2}(map(x->x(t),xr) - θ)
        de = SVector{2}(map(x->x(t),vr) - dθ)
        tau = kp*e + kv*de
        du[1:2] = dθ
        du[3:4] = inv(M)*(tau - C * dθ - G)
        du[5:6] = tau
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(2,sol2)
end

"""
Robô kuka com 7 graus de liberdade na vertical
"""
function kukaRobot(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::Vector{Float64}) where {T<:AbstractMatrix, Y<: AbstractFloat}
    function myrobot(du, u, p, t)
        θ = SVector{7}(u[1:7])
        dθ = SVector{7}(u[8:14])
        set_configuration!(state_kuka,θ)
        set_velocity!(state_kuka,dθ)
        h = SVector{7}(dynamics_bias(state_kuka))
        m = SMatrix{7,7}(mass_matrix(state_kuka))
        e = SVector{7}(xr - θ)
        tau = kp*e - kv*dθ
        tau = SVector{7}(tau)
        du[1:7] = dθ
        du[8:14] = inv(m)*(tau - h)
        du[15:21] = tau
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot, zeros(21), tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(7,sol2)
end

"""
Robô kuka com 7 graus de liberdade na vertical
"""
function kukaRobot(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::Vector{Z}, vr::Vector{Z}) where {T<:AbstractMatrix, Y<: AbstractFloat, Z<:Function}
    function myrobot(du, u, p, t)
        θ = SVector{7}(u[1:7])
        dθ = SVector{7}(u[8:14])
        set_configuration!(state_kuka,θ)
        set_velocity!(state_kuka,dθ)
        h = SVector{7}(dynamics_bias(state_kuka))
        m = SMatrix{7,7}(mass_matrix(state_kuka))
        e = SVector{7}(map(x->x(t),xr) - θ)
        de = SVector{7}(map(x->x(t),vr) - dθ)
        tau = kp*e + kv*de
        tau = SVector{7}(tau)
        du[1:7] = dθ
        du[8:14] = inv(m)*(tau - h)
        du[15:21] = tau
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot, zeros(21), tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(7,sol2)
end


######################## Como utiliza a função ####################
#Ts = 0.08
#tend = 2.0
#t0 = 0.0
#r1 = 1.2
#r2 = 0.6
#kp = SMatrix{2,2}(diagm([10., 30.]))
#kv = SMatrix{2,2}(diagm([5., 3.]))
#x, v, t, a, ta, j, tj = robot2dof(kp, kv, Ts, t0, tend, [sin, cos], [cos, sin])
##################################################################
