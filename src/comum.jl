function organize(mysize::Int, data)
    out_x = map(x -> x[1:mysize],data.u)
    out_dx = map(x -> x[(mysize+1):(2*mysize)],data.u)
    out_int_tau = map(x -> x[((2*mysize)+1):end],data.u)
    const T = diff(data.t)[1]
    #velocidade
    #Aqui estou fazendo uma aproximação da aceleração e do jerk
    #out_d2x = diff(out_dx)/Ts #aceleração
    #out_d3x = diff(out_d2x)/Ts #jerk
    θ = []
    for i = 1:mysize
        push!(θ, map(x -> x[i],out_x))
    end
    ω = []
    for i = 1:mysize
        push!(ω, map(x -> x[i],out_dx))
    end
    ∫τ = []
    for i = 1:mysize
        push!(∫τ, map(x -> x[i],out_int_tau))
    end
    α = []
    for i = 1:mysize
        push!(α, diff(ω[i])./T)
    end
    ta = data.t[1:length(α[1])]
    J = []
    for i = 1:mysize
        push!(J, diff(α[i])./T)
    end
    tj = data.t[1:length(J[1])]
    τ = []
    for i = 1:mysize
        push!(τ,diff(∫τ[i])./T)
    end

    θ, ω, data.t, α, ta, J, tj, τ, ta
end


function tabela(dado::Vector, nome::String)
    size = length(dado)

    # título da tabela
    title = "---"
    for i = 1:size
        title*= "|junta $(i)"
    end
    title*= "\n"

    # barra separadora
    bar = "---"
    for i = 1:size
        bar*= "|---"
    end
    bar *= "\n"

    #valor máximo
    max = "**$(nome) máximo**"
    for i = 1:size
        max*= "| $(round(maximum(abs.(dado[i])),2)) "
    end
    max *= "\n"

    #valor mínimo
    min = "**$(nome) mínimo**"
    for i = 1:size
        min*= "| $(round(minimum(abs.(dado[i])),2)) "
    end
    min *= "\n"

    #valor total
    total = "**$(nome) total**"
    for i = 1:size
        total*= "| $(round(sum(abs.(dado[i])),2)) "
    end
    total *= "\n"
    Markdown.parse(title*bar*max*min*total)
end
tabela([[1.,2.,-3.],[10.,-21.,34.,12.]],"erro")

#Aqui fica a parte de geração de trajetória mínimo jerk

"""
Retorna os coeficientes do polinômio de trajetória com base no princípio do minimum jerk considerando o tempo inicial sendo 0(zero).
Entradas
x0 -> posição inicial
v0 -> velocidade inicial
a0 -> aceleração inicial
t0 -> tempo inicial
xf -> posição final
vf -> velocidade final
af -> aceleração final
tf -> tempo para chegar ao ponto final
"""
function minimumjerk(x0::T, v0::T, a0::T, t0::T, xf::T, vf::T, af::T, tf::T)  where {T<:AbstractFloat}
    b=[ 1 t0 t0^2 t0^3 t0^4 t0^5;
    0 1 2*t0 3*(t0^2) 4*(t0^3) 5*(t0^4);
    0 0 2 6*t0 12*(t0^2) 20*(t0^3);
    1 tf tf^2 tf^3 tf^4 tf^5;
    0 1 2*tf 3*(tf^2) 4*(tf^3) 5*(tf^4);
    0 0 2 6*tf 12*(tf^2) 20*(tf^3)]
    c = [x0, v0, a0, xf, vf, af]
    a = inv(b)*c
    return a
end


"""
Retorna os coeficientes do polinômio de trajetória com base no princípio do minimum jerk considerando o tempo inicial sendo 0(zero).
Entradas
x0 -> posição inicial
v0 -> velocidade inicial
a0 -> aceleração inicial
xf -> posição final
vf -> velocidade final
af -> aceleração final
T -> tempo para chegar ao ponto final
"""
function minimumjerk(x0::T, v0::T, a0::T, xf::T, vf::T, af::T, tf::T) where T<:AbstractFloat
    const t0 = 0
    minimumjerk(x0,v0,a0,t0,xf,vf,af,tf)
end

function minimumjerkf(x0::T, v0::T, a0::T, t0::T, xf::T, vf::T, af::T, tf::T) where T<:AbstractFloat
    functionform(minimumjerk(x0,v0,a0,t0,xf,vf,af,tf))
end

function minimumjerkf(x0::T, v0::T, a0::T, xf::T, vf::T, af::T, tf::T) where T<:AbstractFloat
    const t0 = 0.0
    functionform(minimumjerk(x0,v0,a0,t0,xf,vf,af,tf))
end

# Retorna as funções de posição, velocidade, aceleração e jerk a parti dos coeficientes do polinômio do mínimo jerk
# Entradas
# a -> vetor com os coeficientes do polinômio
function functionform(a::Vector)
    out1 = eval(:(t -> $(a[1]) + $(a[2])*t + $(a[3])*t^2 + $(a[4])*t^3 + $(a[5])*t^4 + $(a[6])*t^5))
    out2 = eval(:(t -> $(a[2]) + 2*$(a[3])*t + 3*$(a[4])*(t^2) + 4*$(a[5])*(t^3) + 5*$(a[6])*(t^4)))
    out3 = eval(:(t ->  2*$(a[3]) + 6*$(a[4])*t + 12*$(a[5])*(t^2) + 20*$(a[6])*(t^3)))
    out4 = eval(:(t -> 6*$(a[4]) + 24*$(a[5])*t + 60*$(a[6])*(t^2)))
    return out1, out2, out3, out4
end


# TODO: fazer a geração de trajetórias pelo método tradicional do livro do craig
