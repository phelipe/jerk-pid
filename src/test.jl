include("comum.jl")
include("modelos.jl")
using Evolutionary, Plots;
pyplot();

Ts     = 0.05 # Intervalo entre leituras da saída
tend   = 2.0  # tempo final para estabilização
t0     = 0.0  # instante inicial
r1     = 0.6  # referência junta 1
r2     = 0.8  # referência junta 2
popul  = 50   # população
iterac = 50;   #iterações

function gerador6(n)
    out = rand(n).*[10000., 100., 1000., 100.]
end

function custo(gain::Vector{Float64})
    kp = SMatrix{2,2}(diagm([gain[1], gain[2]]))
    kv = SMatrix{2,2}(diagm([gain[3], gain[4]]))
    x, v, t, a, ta, j, tj = robot2dof(kp, kv, Ts, t0, tend, [r1, r2])
    erro1 =  - (x[1] - r1)
    erro2 =  - (x[2] - r2)
    sizeVector = length(erro1)

    erro_end_1 = sum(abs.(erro1[floor(Integer,sizeVector/3):end]))
    erro_end_2 = sum(abs.(erro2[floor(Integer,sizeVector/3):end]))

    jerk_1 = sum(abs.(j[1]))
    jerk_2 = sum(abs.(j[2]))

    erro_end = erro_end_1 + erro_end_2
    jerk = jerk_1 + jerk_2

    erro_end = erro_end*10.
    jerk = jerk*0.01

    out =  erro_end + jerk
    #println(" $(erro_end) | $(jerk) | $(out)")
    #println("$(gain)")
    out
end;

N = 4
 result, fitness, cnt = ga(custo, N; initPopulation = gerador6, populationSize = popul, ɛ = 0.1, selection = roulette, crossover = intermediate(0.25), mutation = domainrange(fill(0.5,N)), iterations = iterac)
# result, fitness, cnt = ga(custo, N; initPopulation = gerador6, populationSize = popul, ɛ = 0.1, selection = sus, crossover = intermediate(0.25), mutation = domainrange(fill(0.5,N)), iterations = iterac)
t_end_new = tend
kp = SMatrix{2,2}(diagm(result[1:2]))
kv = SMatrix{2,2}(diagm(result[3:4]))

x, v, t, a, ta, j, tj = robot2dof(kp, kv, Ts, t0, t_end_new, [r1, r2]);

table1 = "|-------- | junta 1  | junta 2 |
|--------| ------------- | ------------- |
|**erro final**| $(rad2deg(x[1][end] - r1)) graus  | $(rad2deg(x[2][end] - r2)) graus  |
|**total jerk **| $(sum(abs.(j[1])))  rad/sec³  | $(sum(abs.(j[2])))  rad/sec³  |
|**máximo jerk**| $(maximum(abs.(j[1])))  rad/sec³  | $(maximum(abs.(j[2])))  rad/sec³ |"

function plotx()
    p1 = plot(t,x[1], label = "PD ótimo - junta 1")
    p1= plot!([r1],seriestype= :hline, label = "referência");
    p2 = plot(t,x[2], label = "PD ótimo - junta 2")
    p2 = plot!([r2],seriestype= :hline, label = "referência");
    plot(p1,p2)
end



function plotj()
    p1 = plot(tj,j[1], label = "jerk 1")
    p2 = plot(tj,j[2], label = "jerk 2")
    plot(p1,p2)
end;


kp_pid = SMatrix{2,2}(diagm([2800., 80.]))
kv_pid = SMatrix{2,2}(diagm([315., 15.]))
x_pid, v_pid, t_pid, a_pid, ta_pid, j_pid, tj_pid = robot2dof(kp_pid, kv_pid, Ts, t0, tend, [r1, r2on
erro1 = -(x_pid[1] - r1)
erro2 = -(x_pid[2] - r2)
table2 = "|-------- | junta 1  | junta 2 |
|--------| ------------- | ------------- |
|**erro final**| $(rad2deg(x_pid[1][end] - r1)) graus  | $(rad2deg(x_pid[2][end] - r2)) graus  |
|**total jerk **| $(sum(abs.(j_pid[1])))  rad/sec³  | $(sum(abs.(j_pid[2])))  rad/sec³  |
|**máximo jerk**| $(maximum(abs.(j_pid[1])))  rad/sec³  | $(maximum(abs.(j_pid[2])))  rad/sec³ |"



#p1 = plot(t_pid,x_pid[1], label = "PD tradicional - junta 1")
#p1= plot!([r1],seriestype= :hline, label = "referência");
#p2 = plot(t_pid,x_pid[2], label = "PD tradicional- junta 2")
#p2 = plot!([r2],seriestype= :hline, label = "referência");
