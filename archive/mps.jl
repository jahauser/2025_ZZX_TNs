# Id = [1 0 
#       0 1]
# PauliX = [0 1
#           1 0]
# PauliY = [0 -1im
#           1im 0]
# PauliZ = [1 0
#           0 -1]

Πx0 = (I+PauliX)/2
Πx1 = (I-PauliX)/2
Πzz0 = (I+kron(PauliZ, PauliZ))/2
Πzz1 = (I-kron(PauliZ, PauliZ))/2
Πxx0 = (I+kron(PauliX, PauliX))/2
Πxx1 = (I-kron(PauliX, PauliX))/2

weak_Πzz0(λ::Float64) = (I + λ*kron(PauliZ,PauliZ))/sqrt(2*(1+λ^2))
weak_Πzz1(λ::Float64) = (I - λ*kron(PauliZ,PauliZ))/sqrt(2*(1+λ^2))

weak_Πx0(λ::Float64) = (I + λ*PauliX)/sqrt(2*(1+λ^2))
weak_Πx1(λ::Float64) = (I - λ*PauliX)/sqrt(2*(1+λ^2))

# function zerostate(L::Int)
#     N = 2L
#     sites = siteinds("Qubit", N)
#     ρ = MPS(sites, _ -> "0")
#     return ρ, sites
# end

# function plusstate(L::Int)
#     N = 2L
#     sites = siteinds("Qubit", N)
#     ρ = MPS(sites, _ -> "+")
#     return ρ, sites
# end

# function ghz(L::Int)
#     N = 2L+2
#     sites = siteinds("Qubit", N)
#     ρ00 = MPS(sites, _ -> "0")
#     ρ11 = MPS(sites, _ -> "1")
#     ρ01 = MPS(sites, i -> mod(i,2)==0 ? "0" : "1")
#     ρ10 = MPS(sites, i -> mod(i,2)==0 ? "1" : "0")
#     return (ρ00 + ρ01 + ρ10 + ρ11)/2, sites
# end

# function x_ghz(L::Int)
#     N = 2L+2
#     sites = siteinds("Qubit", N)
#     ρ00 = MPS(sites, _ -> "+")
#     ρ11 = MPS(sites, _ -> "-")
#     ρ01 = MPS(sites, i -> mod(i,2)==0 ? "+" : "-")
#     ρ10 = MPS(sites, i -> mod(i,2)==0 ? "-" : "+")
#     return (ρ00 + ρ01 + ρ10 + ρ11)/2, sites
# end

function id_mps(s1::Index, s2::Index)
    return MPS([1; 0; 0; 1], [s1, s2])
end

function x_mps(s1::Index, s2::Index)
    return MPS([0; 1; 1; 0], [s1, s2])
end

function z_mps(s1::Index, s2::Index)
    return MPS([1; 0; 0; -1], [s1, s2])
end

function bell(sites::Vector{Index{Int}})
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        a, b = id_mps(sites[i], sites[i+1])
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end

function X_bra(sites::Vector{Index{Int}}, pos::Int)
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        if i == 2pos-1
            a, b = x_mps(sites[i], sites[i+1])
        else
            a, b = id_mps(sites[i], sites[i+1])
        end
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end

function XX_bra(sites::Vector{Index{Int}}, pos1::Int, pos2::Int)
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        if i == 2pos1-1 || i == 2pos2-1
            a, b = x_mps(sites[i], sites[i+1])
        else
            a, b = id_mps(sites[i], sites[i+1])
        end
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end
    
function ZZ_bra(sites::Vector{Index{Int}}, pos1::Int, pos2::Int)
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        if i == 2pos1-1 || i == 2pos2-1
            a, b = z_mps(sites[i], sites[i+1])
        else
            a, b = id_mps(sites[i], sites[i+1])
        end
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end

function ketnorm(ρ::MPS)
    return inner(bell(siteinds(ρ)), ρ)
end

function X_noise_layer(K::Float64, sites::Vector{Index{Int}}; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L
        h = op(PauliX, sites[2pos-1])*op(PauliX, sites[2pos])
        push!(gates, exp(K * h))
    end
    return gates
end

function XX_noise_layer(K::Float64, sites::Vector{Index{Int}}; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L-1
        h = op(PauliX, sites[2pos-1])*op(PauliX, sites[2pos])*op(PauliX, sites[2pos+1])*op(PauliX, sites[2pos+2])
        push!(gates, exp(K * h))
    end
    return gates
end

function ZZ_noise_layer(K::Float64, sites::Vector{Index{Int}}; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L-1
        h = op(PauliZ, sites[2pos-1])*op(PauliZ, sites[2pos])*op(PauliZ, sites[2pos+1])*op(PauliZ, sites[2pos+2])
        push!(gates, exp(K * h))
    end
    return gates
end

function rand_X_strong_measurement_layer(ρ::MPS, sites::Vector{Index{Int}}, p::Float64; reference=false)
    L = length(sites)÷2 - reference
    for pos in 1:L
        if rand() < p
            X_val = inner(X_bra(siteinds(ρ), pos), ρ) / ketnorm(ρ)
            if rand() < (1+X_val)/2
                g1 = op(Πx0, sites[2pos-1])
                g2 = op(Πx0, sites[2pos]) 
            else
                g1 = op(Πx1, sites[2pos-1])
                g2 = op(Πx1, sites[2pos])
            end
            ρ = apply([g1, g2], ρ)
        end
    end
    return ρ
end

function rand_XX_strong_measurement_layer(ρ::MPS, sites::Vector{Index{Int}}, p::Float64; reference=false)
    L = length(sites)÷2 - reference
    for pos in 1:L-1
        if rand() < p
            XX_val = inner(XX_bra(siteinds(ρ), pos, pos+1), ρ) / ketnorm(ρ)
            if rand() < (1+XX_val)/2
                g1 = op(Πxx0, sites[2pos-1], sites[2pos+1])
                g2 = op(Πxx0, sites[2pos], sites[2pos+2])
            else
                g1 = op(Πxx1, sites[2pos-1], sites[2pos+1])
                g2 = op(Πxx1, sites[2pos], sites[2pos+2])
            end
            ρ = apply([g1, g2], ρ)
        end
    end
    return ρ
end

function rand_ZZ_strong_measurement_layer(ρ::MPS, sites::Vector{Index{Int}}, p::Float64; reference=false)
    L = length(sites)÷2 - reference
    for pos in 1:L-1
        if rand() < p
            ZZ_val = inner(ZZ_bra(siteinds(ρ), pos, pos+1), ρ) / ketnorm(ρ)
            if rand() < (1+ZZ_val)/2
                g1 = op(Πzz0, sites[2pos-1], sites[2pos+1])
                g2 = op(Πzz0, sites[2pos], sites[2pos+2])
            else
                g1 = op(Πzz1, sites[2pos-1], sites[2pos+1])
                g2 = op(Πzz1, sites[2pos], sites[2pos+2])
            end
            ρ = apply([g1, g2], ρ)
        end
    end
    return ρ
end

function X_strong_measurement_layer(sites::Vector{Index{Int}}, record::Vector; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L
        if record[pos] == 0
            continue
        elseif record[pos] == 1
            g1 = op(Πx0, sites[2pos-1])
            g2 = op(Πx0, sites[2pos]) 
        elseif record[pos] == -1
            g1 = op(Πx1, sites[2pos-1])
            g2 = op(Πx1, sites[2pos])
        end
        push!(gates, g1)
        push!(gates, g2)
    end
    return gates
end

function XX_strong_measurement_layer(sites::Vector{Index{Int}}, record::Vector; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L-1
        if record[pos] == 0
            continue
        elseif record[pos] == 1
            g1 = op(Πxx0, sites[2pos-1], sites[2pos+1])
            g2 = op(Πxx0, sites[2pos], sites[2pos+2])
        elseif record[pos] == -1
            g1 = op(Πxx1, sites[2pos-1], sites[2pos+1])
            g2 = op(Πxx1, sites[2pos], sites[2pos+2])
        end
        push!(gates, g1)
        push!(gates, g2)
    end
    return gates
end

function ZZ_strong_measurement_layer(sites::Vector{Index{Int}}, record::Vector; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L-1
        if record[pos] == 0
            continue
        elseif record[pos] == 1
            g1 = op(Πzz0, sites[2pos-1], sites[2pos+1])
            g2 = op(Πzz0, sites[2pos], sites[2pos+2])
        elseif record[pos] == -1
            g1 = op(Πzz1, sites[2pos-1], sites[2pos+1])
            g2 = op(Πzz1, sites[2pos], sites[2pos+2])
        end
        push!(gates, g1)
        push!(gates, g2)
    end
    return gates
end

function coherent_info2(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)

    S2_QR = -log2(inner(ρ, ρ))

    a, b = id_mps(sites[N-1], sites[N])
    
    ρ_Q = MPS(ρ[1:N-2])
    ρ_Q[N-2] *= ρ[N-1]*ρ[N]*a*b

    S2_Q = -log2(inner(ρ_Q, ρ_Q))
    return S2_Q - S2_QR
end

function reference_entropy(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)

    A = ITensor(1)
    for i in 1:2:N-2
        B, C = id_mps(sites[i],sites[i+1])
        A = A * B * C * ρ[i] * ρ[i+1]
    end
    ρ_R = A * ρ[N-1] * ρ[N]
    M = Matrix(ρ_R, inds(ρ_R))
    λs = abs.(eigvals(M))
    return sum([-λ*log2(λ) for λ in λs if abs(λ) > 1e-16])
end

function reference_entropy2(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)

    A = ITensor(1)
    for i in 1:2:N-2
        B, C = id_mps(sites[i],sites[i+1])
        A = A * B * C * ρ[i] * ρ[i+1]
    end
    ρ_R = A * ρ[N-1] * ρ[N]
    M = Matrix(ρ_R, inds(ρ_R))
    λs = abs.(eigvals(M))
    return -log2(sum([λ*λ for λ in λs if abs(λ) > 1e-16]))
end

function χZZ(ρ::MPS; reference=false)
    sites = siteinds(ρ)
    N = length(sites)
    L = N÷2 - reference

    ρ = ρ/ketnorm(ρ)

    χ = Float64(L)
    for i in 1:L
        for j in i+1:L
            χ += 2*inner(ZZ_bra(sites, i, j), ρ)^2
        end
    end
    return χ/L
end

function χZZZZ(ρ::MPS; reference=false)
    sites = siteinds(ρ)
    N = length(sites)
    L = N÷2 - reference
    
    ρ = ρ/ketnorm(ρ)

    χ = Float64(L)
    for i in 1:L
        for j in i+1:L
            ZZZZ = op("Z", sites[2i-1]) * op("Z", sites[2i]) * op("Z", sites[2j-1]) * op("Z", sites[2j])
            χ += 2*inner(ρ, apply(ZZZZ, ρ))/inner(ρ, ρ)
        end
    end
    return χ/L
end 

function χXX(ρ::MPS; reference=false)
    sites = siteinds(ρ)
    N = length(sites)
    L = N÷2 - reference

    ρ = ρ/ketnorm(ρ)

    χ = Float64(L)
    for i in 1:L
        for j in i+1:L
            χ += 2*inner(XX_bra(sites, i, j), ρ)^2
        end
    end
    return χ/L
end

function χXXXX(ρ::MPS; reference=false)
    sites = siteinds(ρ)
    N = length(sites)
    L = N÷2 - reference
    
    ρ = ρ/ketnorm(ρ)

    χ = Float64(L)
    for i in 1:L
        for j in i+1:L
            XXXX = op("X", sites[2i-1]) * op("X", sites[2i]) * op("X", sites[2j-1]) * op("X", sites[2j])
            χ += 2*inner(ρ, apply(XXXX, ρ))/inner(ρ, ρ)
        end
    end
    return χ/L
end 



function X_coherent_layer(θ::Float64, sites::Vector{Index{Int}}; reference=false)
    L = length(sites)÷2 - reference
    gates = ITensor[]
    for pos in 1:L
        h1 = op(PauliX, sites[2pos-1])
        h2 = op(PauliX, sites[2pos])
        push!(gates, exp(1im * θ * h1))
        push!(gates, exp(-1im * θ * h2))
    end
    return gates
end

function rand_ZZ_weak_measurement_layer(ρ::MPS, sites::Vector{Index{Int}}, λ::Float64; reference=false)
    L = length(sites)÷2 - reference
    for pos in 1:L-1
        ZZ_val = inner(ZZ_bra(siteinds(ρ), pos, pos+1), ρ) / ketnorm(ρ)
        prob = (1 + 2λ/(1+λ^2)*ZZ_val)/2 
        @assert imag(prob) .< 1e-8
        if rand() < abs(prob)
            g1 = op(weak_Πzz0(λ), sites[mod1(2pos-1,2L)], sites[mod1(2pos+1,2L)])
            g2 = op(weak_Πzz0(λ), sites[mod1(2pos,2L)], sites[mod1(2pos+2,2L)])
        else
            g1 = op(weak_Πzz1(λ), sites[mod1(2pos-1,2L)], sites[mod1(2pos+1,2L)])
            g2 = op(weak_Πzz1(λ), sites[mod1(2pos,2L)], sites[mod1(2pos+2,2L)])
        end
        ρ = apply([g1, g2], ρ)
    end
    return ρ
end

function rand_X_weak_measurement_layer(ρ::MPS, sites::Vector{Index{Int}}, λ::Float64; reference=false)
    L = length(sites)÷2 - reference
    for pos in 1:L
        ZZ_val = inner(ZZ_bra(siteinds(ρ), pos, pos+1), ρ) / ketnorm(ρ)
        prob = (1 + 2λ/(1+λ^2)*ZZ_val)/2 
        @assert imag(prob) .< 1e-8
        if rand() < abs(prob)
            g1 = op(weak_Πx0(λ), sites[2pos-1])
            g2 = op(weak_Πx0(λ), sites[2pos])
        else
            g1 = op(weak_Πx1(λ), sites[2pos-1])
            g2 = op(weak_Πx1(λ), sites[2pos])
        end
        ρ = apply([g1, g2], ρ)
    end
    return ρ
end