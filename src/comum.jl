function organize(mysize::Int, data)
    out_x = map(x -> x[1:mysize],data.u)
    out_dx = map(x -> x[(mysize+1):end],data.u)
    const T = diff(data.t)[1]
    #velocidade
    #Aqui estou fazendo uma aproximação da aceleração e do jerk
    #out_d2x = diff(out_dx)/Ts #aceleração
    #out_d3x = diff(out_d2x)/Ts #jerk
    θ = []
    for i=1:mysize
        push!(θ, map(x -> x[i],out_x))
    end
    ω = []
    for i=1:mysize
        push!(ω, map(x -> x[i],out_dx))
    end
    α = []
    for i=1:mysize
        push!(α, diff(ω[i])./T)
    end
    ta = data.t[1:length(α[1])]
    J = []
    for i=1:mysize
        push!(J, diff(α[i])./T)
    end
    tj = data.t[1:length(J[1])]

    θ, ω, data.t, α, ta, J, tj
end
