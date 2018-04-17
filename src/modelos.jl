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
function robot2dof(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::AbstractVector) where {T<:AbstractMatrix, Y<: AbstractFloat}

    function myrobot(du, u, p, t)
        m = SVector{2}([23.902, 1.285])
        l= SVector{2}([.45, .45])
        r= SVector{2}([0.091, 0.048])
        I= SVector{2}([1.266, 0.093])
        g = 9.81
        θ = SVector{2}(u[1:2])
        dθ = SVector{2}(u[3:4])
        const cos1 = cos(θ[1])
        const cos2 = cos(θ[2])
        const sin1 = sin(θ[1])
        const sin2 = sin(θ[2])
        const cos12 = cos(θ[1]+θ[2])

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

        e = xr - θ
        tau = kp*e - kv*dθ
        du[1:2] = dθ
        du[3:4] = inv(M)*(tau - C * dθ - G)
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(2,sol2)
end

"""
Robô kuka com 7 graus de liberdade na vertical
"""
function kukaRobot(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::AbstractVector) where {T<:AbstractMatrix, Y<: AbstractFloat}
    function myrobot(du, u, p, t)
        θ = SVector{7}(u[1:7])
        dθ = SVector{7}(u[8:end])
        set_configuration!(state_kuka,θ)
        set_velocity!(state_kuka,dθ)
        h = SVector{7}(dynamics_bias(state_kuka))
        m = SMatrix{7,7}(mass_matrix(state_kuka))
        e = xr - θ
        tau = kp*e - kv*dθ
        du[1:7] = dθ
        du[8:end] = inv(m)*(tau - h)
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    #,reltol=1e-2,abstol=1e-3,force_dtmin=true
    organize(7,sol2)
end


"""
Robô de estrutura paralela cadeia fechada na horizontal
"""
function parallelRobot(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, xr::AbstractVector) where {T<:AbstractMatrix, Y<: AbstractFloat}

    function myrobot(du, u, p, t)
        #Entradas ativas
        q11 = u[1]
        q21 = u[2]
        dq11 = u[3]
        dq21 = u[4]

        #####################POSIÇÕES PASSIVAS#################
        #comprimento dos links em metros
        a11 = -68/2000
        a12 = 100/1000
        a13 = 125/1000
        a21 = 68/2000
        a22 = 80/1000
        a23 = 115/1000

        C = (a21-a11)+a22*cos(q21)-a12*cos(q11)
        D = a22*sin(q21)-a12*sin(q11)
        E = (a13^2-a23^2-(C^2+D^2))/(2*a23)

        q22 = acos(2*E*C-sqrt((2*E*C)^2-4*(C^2+D^2)*(E^2-D^2))/2*(C^2+D^2))-q21
        q12 = asin(a22*sin(q21)+a23*sin(q21+q22)-a12*sin(q11)/a13)-q11


        ####################VELOCIDADES PASSIVAS#################
        ro11 = 1
        ro12 = 0
        ro21 = 0
        ro22 = 1
        ro31 = -(8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/sin(q11 + q12 - q21 - q22)
        ro32 = -(16*sin(q22))/(25*sin(q11 + q12 - q21 - q22))
        ro41 = (20*sin(q12))/(23*sin(q11 + q12 - q21 - q22))
        ro42 = -(200*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22))

        ro = SMatrix{4,2}([ro11 ro12;ro21 ro22;ro31 ro32;ro41 ro42])

        dq11 ,dq21, dq12, dq22 = ro * SVector{2}([dq11,dq21])


        ############ Dinâmica#############
        D11 = (103109*cos(q12))/200000000 + (217704077823*sin(q12)^2)/(82656250000000*sin(q11 + q12 - q21 - q22)^2) - ((sin(q11 + q12 - q21 - q22) + (4*sin(q11 - q21 - q22))/5)*((103109*cos(q12))/400000000 - (2479099*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/(1500000000*sin(q11 + q12 - q21 - q22)) + 2371507/16000000000))/sin(q11 + q12 - q21 - q22) - (8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10)*((103109*cos(q12))/400000000 + 2371507/16000000000))/sin(q11 + q12 - q21 - q22) + 2177040973243/3000000000000000

        D12 = (20*sin(q12)*((176952321*cos(q22))/62500000000 + 756825076917/250000000000000))/(23*sin(q11 + q12 - q21 - q22)) - (16*sin(q22)*((103109*cos(q12))/400000000 - (2479099*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/(1500000000*sin(q11 + q12 - q21 - q22)) + 2371507/16000000000))/(25*sin(q11 + q12 - q21 - q22)) - (217704077823*sin(q12)*(23*sin(q11 + q12 - q21 - q22) + 16*sin(q11 + q12 - q21)))/(1653125000000000*sin(q11 + q12 - q21 - q22)^2)

        D21 = (20*sin(q12)*((176952321*cos(q22))/62500000000 - (217704077823*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(7187500000000*sin(q11 + q12 - q21 - q22)) + 756825076917/250000000000000))/(23*sin(q11 + q12 - q21 - q22)) - (16*sin(q22)*((103109*cos(q12))/400000000 + 2371507/16000000000))/(25*sin(q11 + q12 - q21 - q22)) + (2479099*sin(q22)*(sin(q11 + q12 - q21 - q22) + (4*sin(q11 - q21 - q22))/5))/(18750000000*sin(q11 + q12 - q21 - q22)^2)

        D22 = (176952321*cos(q22))/31250000000 + (2479099*sin(q22)^2)/(29296875000*sin(q11 + q12 - q21 - q22)^2) - (200*((176952321*cos(q22))/62500000000 + 756825076917/250000000000000)*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22)) - ((23*sin(q11 + q12 - q21 - q22) + 16*sin(q11 + q12 - q21))*((176952321*cos(q22))/62500000000 - (217704077823*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(7187500000000*sin(q11 + q12 - q21 - q22)) + 756825076917/250000000000000))/(23*sin(q11 + q12 - q21 - q22)) + 89825638637/15625000000000

        C11 = (103109*dq12*sin(q12)*(sin(q11 + q12 - q21 - q22) + (4*sin(q11 - q21 - q22))/5))/(400000000*sin(q11 + q12 - q21 - q22)) - ((103109*cos(q12))/400000000 - (2479099*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/(1500000000*sin(q11 + q12 - q21 - q22)) + 2371507/16000000000)*((8*sin(q21 + q22)*((sin(q11 + q12)*(dq11 + dq12))/8 + (dq11*sin(q11))/10))/sin(q11 + q12 - q21 - q22) + (8*cos(q21 + q22)*((cos(q11 + q12)*(dq11 + dq12))/8 + (dq11*cos(q11))/10))/sin(q11 + q12 - q21 - q22) - (20*sin(q12)*((23*cos(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)) + (23*sin(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22))))/(23*sin(q11 + q12 - q21 - q22)) - (8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10)*((cos(q11 + q12)*cos(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22)))/sin(q11 + q12 - q21 - q22)) - (217704077823*sin(q12)*((200*sin(q11 + q12)*((sin(q11 + q12)*(dq11 + dq12))/8 + (dq11*sin(q11))/10))/(23*sin(q11 + q12 - q21 - q22)) - (8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10)*((25*cos(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22)) + (25*sin(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22))))/sin(q11 + q12 - q21 - q22) - (20*sin(q12)*((cos(q11 + q12)*cos(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22)))/(23*sin(q11 + q12 - q21 - q22)) + (200*cos(q11 + q12)*((cos(q11 + q12)*(dq11 + dq12))/8 + (dq11*cos(q11))/10))/(23*sin(q11 + q12 - q21 - q22))))/(71875000000000*sin(q11 + q12 - q21 - q22)) - (103109*dq12*sin(q12))/200000000 - (103109*dq11*sin(q12)*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/(50000000*sin(q11 + q12 - q21 - q22))

        C12 = ((103109*cos(q12))/400000000 - (2479099*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10))/(1500000000*sin(q11 + q12 - q21 - q22)) + 2371507/16000000000)*((8*sin(q21 + q22)*((23*sin(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*sin(q21))/25))/sin(q11 + q12 - q21 - q22) - (200*((23*cos(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)) + (23*sin(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)))*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (16*sin(q22)*((cos(q11 + q12)*cos(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22)))/(25*sin(q11 + q12 - q21 - q22)) + (8*cos(q21 + q22)*((23*cos(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*cos(q21))/25))/sin(q11 + q12 - q21 - q22)) + (217704077823*sin(q12)*((200*sin(q11 + q12)*((23*sin(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*sin(q21))/25))/(23*sin(q11 + q12 - q21 - q22)) - (200*((cos(q11 + q12)*cos(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22))*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (200*cos(q11 + q12)*((23*cos(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*cos(q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (16*sin(q22)*((25*cos(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22)) + (25*sin(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22))))/(25*sin(q11 + q12 - q21 - q22))))/(71875000000000*sin(q11 + q12 - q21 - q22)) + (103109*dq12*sin(q12)*sin(q22))/(625000000*sin(q11 + q12 - q21 - q22)) + (176952321*dq21*sin(q12)*sin(q22))/(71875000000*sin(q11 + q12 - q21 - q22))

        C21 = (2479099*sin(q22)*((8*sin(q21 + q22)*((sin(q11 + q12)*(dq11 + dq12))/8 + (dq11*sin(q11))/10))/sin(q11 + q12 - q21 - q22) + (8*cos(q21 + q22)*((cos(q11 + q12)*(dq11 + dq12))/8 + (dq11*cos(q11))/10))/sin(q11 + q12 - q21 - q22) - (20*sin(q12)*((23*cos(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)) + (23*sin(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22))))/(23*sin(q11 + q12 - q21 - q22)) - (8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10)*((cos(q11 + q12)*cos(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22)))/sin(q11 + q12 - q21 - q22)))/(18750000000*sin(q11 + q12 - q21 - q22)) - ((176952321*cos(q22))/62500000000 - (217704077823*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(7187500000000*sin(q11 + q12 - q21 - q22)) + 756825076917/250000000000000)*((200*sin(q11 + q12)*((sin(q11 + q12)*(dq11 + dq12))/8 + (dq11*sin(q11))/10))/(23*sin(q11 + q12 - q21 - q22)) - (8*(sin(q11 + q12 - q21 - q22)/8 + sin(q11 - q21 - q22)/10)*((25*cos(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22)) + (25*sin(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22))))/sin(q11 + q12 - q21 - q22) - (20*sin(q12)*((cos(q11 + q12)*cos(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22)))/(23*sin(q11 + q12 - q21 - q22)) + (200*cos(q11 + q12)*((cos(q11 + q12)*(dq11 + dq12))/8 + (dq11*cos(q11))/10))/(23*sin(q11 + q12 - q21 - q22))) - (103109*dq11*sin(q12)*sin(q22))/(625000000*sin(q11 + q12 - q21 - q22)) - (176952321*dq22*sin(q12)*sin(q22))/(71875000000*sin(q11 + q12 - q21 - q22))

        C22 = ((176952321*cos(q22))/62500000000 - (217704077823*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(7187500000000*sin(q11 + q12 - q21 - q22)) + 756825076917/250000000000000)*((200*sin(q11 + q12)*((23*sin(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*sin(q21))/25))/(23*sin(q11 + q12 - q21 - q22)) - (200*((cos(q11 + q12)*cos(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq21 + dq22))/sin(q11 + q12 - q21 - q22))*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (200*cos(q11 + q12)*((23*cos(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*cos(q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (16*sin(q22)*((25*cos(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22)) + (25*sin(q11 + q12)^2*(dq11 + dq12))/(23*sin(q11 + q12 - q21 - q22))))/(25*sin(q11 + q12 - q21 - q22))) - (176952321*dq22*sin(q22))/31250000000 - (2479099*sin(q22)*((8*sin(q21 + q22)*((23*sin(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*sin(q21))/25))/sin(q11 + q12 - q21 - q22) - (200*((23*cos(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)) + (23*sin(q21 + q22)^2*(dq21 + dq22))/(25*sin(q11 + q12 - q21 - q22)))*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(23*sin(q11 + q12 - q21 - q22)) + (16*sin(q22)*((cos(q11 + q12)*cos(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22) + (sin(q11 + q12)*sin(q21 + q22)*(dq11 + dq12))/sin(q11 + q12 - q21 - q22)))/(25*sin(q11 + q12 - q21 - q22)) + (8*cos(q21 + q22)*((23*cos(q21 + q22)*(dq21 + dq22))/200 + (2*dq21*cos(q21))/25))/sin(q11 + q12 - q21 - q22)))/(18750000000*sin(q11 + q12 - q21 - q22)) + (176952321*dq22*sin(q22)*(23*sin(q11 + q12 - q21 - q22) + 16*sin(q11 + q12 - q21)))/(1437500000000*sin(q11 + q12 - q21 - q22)) - (176952321*dq21*sin(q22)*((23*sin(q11 + q12 - q21 - q22))/200 + (2*sin(q11 + q12 - q21))/25))/(7187500000*sin(q11 + q12 - q21 - q22))

        D = SMatrix{2,2}([D11 D12;D21 D22])
        C = SMatrix{2,2}([C11 C12;C21 C22])


        θ = SVector{2}([q11, q21])

        dθ = SVector{2}([dq11, dq21])
        e = xr - θ
        tau = kp*e - kv*dθ
        du[1:2] = dθ
        du[3:4] = inv(D)*(tau - C * dθ)

    end
    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts, maxiters = 1e7, force_dtmin=true,reltol=1e-2,abstol=1e-3)
    organize(2,sol2)
end



######################## Como utiliza a função ####################
#Ts = 0.08
#tend = 2.0
#t0 = 0.0
#r1 = 1.2
#r2 = 0.6
#kp = SMatrix{2,2}(diagm([10., 30.]))
#kv = SMatrix{2,2}(diagm([5., 3.]))
#x, v, t, a, ta, j, tj = robot2dof(kp, kv, Ts, t0, tend, r1, r2);
##################################################################
