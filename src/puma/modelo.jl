using StaticArrays
using DifferentialEquations
using RigidBodyDynamics

urdf_puma = "puma560_description/puma560.urdf"
mecanismo_puma = parse_urdf(Float64, urdf_puma)
state_puma = MechanismState(mecanismo_puma)

"""
Robô PUMA com 5 graus de liberdade(contando somente as juntas rotacionais) na vertical
"""
function pumaRobot(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::AbstractVector) where {T<:AbstractMatrix, Y<: AbstractFloat}
    function myrobot(du, u, p, t)
        θ = SVector{5}(u[1:5])
        dθ = SVector{5}(u[6:10])
        set_configuration!(state_puma,θ)
        set_velocity!(state_puma,dθ)
        h = SVector{5}(dynamics_bias(state_puma))
        m = SMatrix{5,5}(mass_matrix(state_puma))
        e = xr - θ
        tau = kp*e - kv*dθ
        tau = SVector{5}(tau)
        du[1:5] = dθ
        du[6:10] = inv(m)*(tau - h)
        du[11:15] = tau
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot, zeros(15), tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(5,sol2)
end
