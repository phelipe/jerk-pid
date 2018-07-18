include("../comum.jl")
include("../modelos.jl")
using Plots
Plots.scalefontsizes(1.2)
pyplot()
#fnt = Plots.font("Helvetica", 10.0)
#default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

Ts     = 0.05 # Intervalo entre leituras da saída
tend   = 2.0  # tempo final para estabilização
t0     = 0.0  # instante inicial

popul  = 20   # população
iterac = 15   #iterações
α = 10.       #parâmetro para o erro
β = 0.01      #parâmetro para o jerk
γ = 0.1       #parâmetro para o torque
per = 1/2     #inicio da leitura do vetor a parti de per do comprimento total

x1, v1, a1, j1 = minimumjerkf(0.0, 0.0, 0.0,0.0, 0.6, 0.0, 0.0, tend)
x2, v2, a2, j2 = minimumjerkf(0.0, 0.0, 0.0,0.0, 0.8, 0.0, 0.0, tend)
xr = [x1,x2]
vr = [v1,v2]
ar = [a1, a2]
jr = [j1,j2];

#Dados otimizados
ganhos = [6748.16, 2996.68, 969.85, 66.5882]
kp = SMatrix{2,2}(diagm(ganhos[1:2]))
kv = SMatrix{2,2}(diagm(ganhos[3:4]))
x, v, t, a, ta, j, tj, τ, t_tau = robot2dof(kp, kv, Ts, t0, tend, xr, vr)
erro1 = -(x[1] .- map(i->xr[1](i), t))
erro2 = -(x[2] .- map(i->xr[2](i), t))
erro1a = -(a[1] .- map(i->ar[1](i), ta))
erro2a = -(a[2] .- map(i->ar[2](i), ta))
erro1j = -(j[1] .- map(i->jr[1](i), tj))
erro2j = -(j[2] .- map(i->jr[2](i), tj))
erro = [erro1, erro2]

#Dados PD
kp_pid = SMatrix{2,2}(diagm([8550., 160.]))
kv_pid = SMatrix{2,2}(diagm([415., 15.]))
x_pid, v_pid, t_pid, a_pid, ta_pid, j_pid, tj_pid, τ_pid, t_tau_pid = robot2dof(kp_pid, kv_pid, Ts, t0, tend, xr,vr)
erro1_pd = -(x_pid[1] .- map(i->xr[1](i), t_pid))
erro2_pd = -(x_pid[2] .- map(i->xr[2](i), t_pid))
erro1a_pd = -(a_pid[1] .- map(i->ar[1](i), ta_pid))
erro2a_pd = -(a_pid[2] .- map(i->ar[2](i), ta_pid))
erro1j_pd = -(j_pid[1] .- map(i->jr[1](i), tj_pid))
erro2j_pd = -(j_pid[2] .- map(i->jr[2](i), tj_pid))

function plotx(aux)
    p1 = plot(t,x[aux], label = " MJ-OPD  - joint $(aux)",xlabel ="Time (s)", ylabel = "Position (rad)", line=(2,:dot))
    p1 = plot!(t_pid,x_pid[aux], label ="ZN-PD - joint $(aux)",line=(2,:dash))
    p1= plot!(t,map(i->xr[aux](i), t), label = "Desired",line=(1))
    plot(p1)
end

function plotj(aux)
    p1 = plot(tj,j[aux], label = " MJ-OPD  - joint $(aux)", xlabel ="Time (s)", ylabel = "Jerk (rad/s³)",line=(2,:dot))
    p1 = plot!(tj_pid,j_pid[aux], label ="ZN-PD - joint $(aux)",line=(2,:dash))
    p1= plot!(t,map(i->jr[aux](i), t), label = "Desired",line=(1))
    plot(p1)
end;

function plotTau(aux)
    p1 = plot(t_tau,τ[aux], label = " MJ-OPD  - joint $(aux)", xlabel ="Time (s)", ylabel = "Torque (Nm)",line=(2,:dot))
    p1 = plot!(t_tau_pid,τ_pid[aux], label ="ZN-PD - joint $(aux)",line=(2,:dash))
    plot(p1)
end;

function plotv(aux)
    p1 = plot(t,v[aux], label = "  MJ-OPD - joint $(aux)", xlabel ="Time (s)", ylabel = "Velocity (rad/s)",line=(2,:dot))
    p1 = plot!(t_pid,v_pid[aux], label ="ZN-PD - joint $(aux)",line=(2,:dash))
    p1= plot!(t,map(i->vr[aux](i), t), label = "Desired",line=(1))
    plot(p1)
end;

function plota(aux)
    p1 = plot(ta,a[aux], label = " MJ-OPD  - joint $(aux)", xlabel ="Time (s)", ylabel = "Acceleration (rad/s²)",line=(2,:dot))
    p1 = plot!(ta_pid,a_pid[aux], label ="ZN-PD - joint $(aux)",line=(2,:dash))
    p1= plot!(t,map(i->ar[aux](i), t), label = "Desired",line=(1))
    plot(p1)
end;
