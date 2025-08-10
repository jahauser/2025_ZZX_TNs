# works for non-normalized
function doubled_expval(ρ::MPS, M::AbstractMatrix, pos::Int; ref=false)
    sites = siteinds(ρ)

    bra = M_bra(sites, M, pos; ref=ref)
    val = inner(bra, ρ) / doubledtrace(ρ)
    return val
end

function doubled_measure(ρ::MPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; ref=false)
    sites = siteinds(ρ)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - ref

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g1 = op(Π, [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]...)
    g2 = op(Π, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...)

    ρ = apply([g1, g2], ρ)
    return ρ
end

function doubled_measure(ρ::MPS, M::AbstractMatrix, λ::Float64, pos::Int; ref=false)
    val = doubled_expval(ρ, M, pos; ref=ref)
    prob = (1 + 2λ/(1+λ^2)*val)/2  

    if rand() < abs(prob)
        return doubled_measure(ρ, M, λ, pos, false; ref=ref), false, val
    else
        return doubled_measure(ρ, M, λ, pos, true; ref=ref), true, val
    end
end

function forced_doubled_measure_with_val(ρ::MPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; ref=false)
    return doubled_measure(ρ, M, λ, pos, m; ref=ref), m, doubled_expval(ρ, M, pos; ref=ref)
end

function doubled_measure(ρ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; ref=false)
    for (i,pos) in enumerate(positions)
        ρ = doubled_measure(ρ, M, λ, pos, ms[i]; ref=ref)
    end
    return ρ
end

function doubled_measure(ρ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector; ref=false)
    ms = zeros(Bool, length(positions))
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        ρ, m, val = doubled_measure(ρ, M, λ, pos; ref=ref)
        ms[i] = m
        vals[i] = val
    end
    return ρ, ms, vals
end

function forced_doubled_measure_with_val(ρ::MPS, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; ref=false)
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        ρ, _, val = forced_doubled_measure_with_val(ρ, M, λ, pos, ms[i]; ref=ref)
        vals[i] = val
    end
    return ρ, ms, vals
end