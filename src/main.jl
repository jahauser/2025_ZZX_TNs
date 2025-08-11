Id = [1 0 
      0 1]

PauliX = [0 1
          1 0]

PauliY = [0 -1im
          1im 0]

PauliZ = [1 0
          0 -1]

SWAP = [1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1]


ITensors.op(::OpName"ZZ",::SiteType"Doubled") = [1 0 0 0
                                                 0 -1 0 0
                                                 0 0 -1 0
                                                 0 0 0 1]


function combined(ρ::MPS)
    sites = siteinds(ρ)
    N = length(sites)
    return MPS([combiner(sites[i:i+1]; tags="Doubled,$i") * (ρ[i] * ρ[i+1]) for i in 1:2:N])
end

# TODO: could someday adjust to take out exponential so p = 0.5 works better
# assumes M satisfies M^2 = 1
# function decoherence_layer(sites::Vector{Index{Int}}, M::AbstractMatrix, p::Float64, positions::AbstractVector, L::Int)
#     K = atanh(p/(1-p))
#     M_width = Int(log2(size(M)[1]))

#     gates = ITensor[]
#     for pos in positions
#         # println([mod1(2*(pos+i)-1,2L) for i in 0:M_width-1])
#         # println([mod1(2*(pos+i),2L) for i in 0:M_width-1])
#         legs1 = [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]
#         legs2 = [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]
#         h = op(M, legs1...)*op(M, legs2...)
#         push!(gates, exp(K * h))
#     end
#     return gates
# end

function decoherence_layer(sites::Vector{Index{Int}}, M::AbstractMatrix, p::Float64, positions::AbstractVector, L::Int)
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        # println([mod1(2*(pos+i)-1,2L) for i in 0:M_width-1])
        # println([mod1(2*(pos+i),2L) for i in 0:M_width-1])
        legs1 = [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]
        legs2 = [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]
        g = p*op(I, legs1...)*op(I, legs2...) + (1-p)*op(M, legs1...)*op(M, legs2...)
        push!(gates, g)
    end
    return gates
end

function coherent_layer(sites::Vector{Index{Int}}, M::AbstractMatrix, θ::Float64, positions::AbstractVector, L::Int)
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        # println([mod1(2*(pos+i)-1,2L) for i in 0:M_width-1])
        # println([mod1(2*(pos+i),2L) for i in 0:M_width-1])
        legs1 = [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]
        legs2 = [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]
        h1 = op(exp(1im*θ*M), legs1...)
        h2 = op(exp(-1im*θ*M), legs2...)
        push!(gates, h1)
        push!(gates, h2)
    end
    return gates
end

function singled_coherent_layer(sites::Vector{Index{Int}}, M::AbstractMatrix, θ::Float64, positions::AbstractVector, L::Int)
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        legs = [sites[mod1(pos+i,L)] for i in 0:M_width-1]
        h = op(exp(1im*θ*M), legs...)
        push!(gates, h)
    end
    return gates
end

function id_mps(s1::Index, s2::Index)
    return MPS([1; 0; 0; 1], [s1, s2])
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

function M_bra(sites::Vector{Index{Int}}, M::AbstractMatrix, pos::Int; ref=false)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - ref

    bra = bell(sites)
    bra = apply(op(M, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...), bra)
    return bra
end

function doubledtrace(ρ::MPS)
    return inner(bell(siteinds(ρ)), ρ)
end

function Svn(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    _, Λ, _ = svd(ψ[b], linkind(ψ, b))

    ps = [Λ[i,i]^2 for i in 1:dim(Λ,1) if Λ[i,i] > 0]
    return -sum([p*log(p) for p in ps])
end