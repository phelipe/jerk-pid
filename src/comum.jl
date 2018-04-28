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
