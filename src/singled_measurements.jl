# works for non-normalized
function singled_expval(ψ::MPS, M::AbstractMatrix, pos::Int; ref=false)
    sites = siteinds(ψ)
    L = length(sites) - ref
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ)
    val = inner(ψ, Mψ) / norm(ψ)^2
    return val
end

function singled_measure(ψ::MPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; ref=false)
    sites = siteinds(ψ)
    L = length(sites) - ref
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g = op(Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ = apply(g, ψ)
    return ψ
end

function singled_measure(ψ::MPS, M::AbstractMatrix, λ::Float64, pos::Int; ref=false)
    val = singled_expval(ψ, M, pos; ref=ref)
    prob = (1 + 2λ/(1+λ^2)*val)/2  

    if rand() < abs(prob)
        return singled_measure(ψ, M, λ, pos, false; ref=ref), false, val
    else
        return singled_measure(ψ, M, λ, pos, true; ref=ref), true, val
    end
end

function forced_singled_measure_with_val(ψ::MPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; ref=false)
    return singled_measure(ψ, M, λ, pos, m; ref=ref), m, singled_expval(ψ, M, pos; ref=ref)
end

function singled_measure(ψ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; ref=false)
    for (i,pos) in enumerate(positions)
        ψ = singled_measure(ψ, M, λ, pos, ms[i]; ref=ref)
    end
    return ψ
end

function singled_measure(ψ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector; ref=false)
    ms = zeros(Bool, length(positions))
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        ψ, m, val = singled_measure(ψ, M, λ, pos; ref=ref)
        ms[i] = m
        vals[i] = val
    end
    return ψ, ms, vals
end

function forced_singled_measure_with_val(ψ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; ref=false)
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        ψ, _, val = forced_singled_measure_with_val(ψ, M, λ, pos, ms[i]; ref=ref)
        vals[i] = val
    end
    return ψ, ms, vals
end

