function ghz(L::Int; ref=false)
    N = 2L+2ref
    sites = siteinds("Qubit", N)
    ρ00 = MPS(sites, _ -> "0")
    ρ11 = MPS(sites, _ -> "1")
    ρ01 = MPS(sites, i -> mod(i,2)==0 ? "0" : "1")
    ρ10 = MPS(sites, i -> mod(i,2)==0 ? "1" : "0")
    return (ρ00 + ρ01 + ρ10 + ρ11)/2, sites
end

function pure_ghz(L::Int; ref=false)
    N = L+ref
    sites = siteinds("Qubit", N)
    ψ0 = MPS(sites, _ -> "0")
    ψ1 = MPS(sites, _ -> "1")
    return (ψ0 + ψ1)/sqrt(2), sites
end


function to_matrix(ρ::MPS)
    ρ = contract(ρ)
    tensor_inds = inds(ρ)
    permutation = vcat(reverse(tensor_inds[2:2:end]), reverse(tensor_inds[1:2:end]))

    M = Array{ComplexF64,length(tensor_inds)}(ρ, permutation)

    return reshape(M, (2^(length(tensor_inds)÷2), 2^(length(tensor_inds)÷2)))
end