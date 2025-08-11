# focused on open boundary conditions for now
function ED_ghzstate(L::Int; ref=false)
    ψ = spzeros(2^(L+ref))
    ψ[1] = 1/sqrt(2)
    ψ[end] = 1/sqrt(2)
    return ψ
end

function densitymatrix(ψ::AbstractVector)
    return ψ * ψ'
end

dm = densitymatrix

ED_ghzdm = dm ∘ ED_ghzstate

function pad(M::AbstractMatrix, before::Int, after::Int)
    return kron(sparse(I, 2^before, 2^before), M, sparse(I, 2^after, 2^after))
end

# assumes M satisfies M^2 = 1
function decohere(ρ::AbstractMatrix, M::AbstractMatrix, p::Float64, site::Int)
    L = Int(log2(size(ρ)[1]))
    M_width = Int(log2(size(M)[1]))
    padded_M = pad(M, site-1, L+1-site-M_width)
    return (1-p)*ρ + p*padded_M*ρ*padded_M
end

function decohere(ρ::AbstractMatrix, M::AbstractMatrix, p::Float64, sites::AbstractVector)
    for site in sites
        ρ = decohere(ρ, M, p, site)
    end
    return ρ
end

function expval(ψ::AbstractVector, M::AbstractMatrix, site::Int; classical=false)
    L = Int(log2(size(ψ)[1]))
    M_width = Int(log2(size(M)[1]))
    if classical
        padded_M = kron(ones(2^(site-1)), M*ones(2^M_width), ones(2^(L+1-site-M_width)))
        return (padded_M' * ψ)/(sum(ψ))
    else
        padded_M = pad(M, site-1, L+1-site-M_width)
        return (ψ' * padded_M * ψ)/norm(ψ)^2
    end
end

function expval(ρ::AbstractMatrix, M::AbstractMatrix, site::Int)
    L = Int(log2(size(ρ)[1]))
    M_width = Int(log2(size(M)[1]))
    padded_M = pad(M, site-1, L+1-site-M_width)
    return tr(ρ * padded_M)/tr(ρ)
end

# assumes M satisfies M^2 = 1
function measure(ρ::AbstractMatrix, M::AbstractMatrix, λ::Float64, site::Int, m::Bool)
    L = Int(log2(size(ρ)[1]))
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))

    padded_Π = pad(Π, site-1, L+1-site-M_width)
    ρ =  padded_Π*ρ*padded_Π
    return ρ/tr(ρ)
end

function measure(ψ::AbstractVector, M::AbstractMatrix, λ::Float64, site::Int, m::Bool)
    L = Int(log2(size(ψ)[1]))
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))

    padded_Π = pad(Π, site-1, L+1-site-M_width)
    ψ =  padded_Π*ψ
    return ψ/sqrt(ψ'*ψ)
end

function measure(ϕ::AbstractArray, M::AbstractMatrix, λ::Float64, site::Int; classical=false)
    val = expval(ϕ, M, site; classical=classical)
    prob = (1 + 2λ/(1+λ^2)*val)/2
    if rand() < prob
        return measure(ϕ, M, λ, site, false), false, val
    else
        return measure(ϕ, M, λ, site, true), true, val
    end
end

function measure(ϕ::AbstractArray, M::AbstractMatrix, λ::Float64, sites::AbstractVector, ms::Vector{Bool})
    for (i,site) in enumerate(sites)
        ϕ = measure(ϕ, M, λ, site, ms[i])
    end 
    return ϕ
end

function forced_measure_with_val(ϕ::AbstractArray, M::AbstractMatrix, λ::Float64, site::Int, m::Bool)
    return measure(ϕ, M, λ, site, m), m, expval(ϕ, M, site)
end

function forced_measure_with_val(ϕ::AbstractArray, M::AbstractMatrix, λ::Float64, sites::AbstractVector, ms::Vector{Bool})
    vals = zeros(Float64, length(sites))
    for (i,site) in enumerate(sites)
        ϕ, _, val = forced_measure_with_val(ϕ, M, λ, site, ms[i])
        vals[i] = val
    end 
    return ϕ, ms, vals
end

function measure(ϕ::AbstractArray, M::AbstractMatrix, λ::Float64, sites::AbstractVector{Int}; classical=false)
    ms = zeros(Bool, length(sites))
    vals = zeros(Float64, length(sites))

    for (i,site) in enumerate(sites)
        ϕ, m, val = measure(ϕ, M, λ, site; classical=classical)
        ms[i] = m
        vals[i] = val
    end 
    return ϕ, ms, vals
end

function ED_full(L::Int, T::Int, λ::Float64, δ::Float64, q::Float64, X_ms::Matrix{Bool}, ZZ_ms::Matrix{Bool}; observables=Vector{Symbol}())
    ρ = ED_ghzdm(L; ref=true)

    data = Dict([s => zeros(Float64, T+1) for s in observables])
    # data = ED_update_data(ρ, observables, data, 1)
    ZZ_vals = zeros(Float64, T, L-1)
    X_vals = zeros(Float64, T, L)

    λzz = δ*(1-λ)
    λx = δ*λ

    for t in 1:T
        ρ = decohere(ρ, PauliX, q, 1:L)
        
        ρ, _, X_vals[t,:] = forced_measure_with_val(ρ, PauliX, λx, 1:L, X_ms[t,:])

        ρ = decohere(ρ, kron(PauliZ,PauliZ), q, 1:L-1)
        
        ρ, _, ZZ_vals[t,:] = forced_measure_with_val(ρ, kron(PauliZ,PauliZ), λzz, 1:L-1, ZZ_ms[t,:])

        # data = ED_update_data(ρ, observables, data, t+1)
    end
    return ρ, (X_vals, ZZ_vals), data
end