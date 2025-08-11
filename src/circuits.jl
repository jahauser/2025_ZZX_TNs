function MPS_update_data(ρ::MPS, observables::Vector{Symbol}, data::Dict{Symbol,Vector{Float64}}, t::Int; ref=false)
    if :Ic in observables
        data[:Ic][t] = Ic2(ρ)
    end
    if :SR in observables
        data[:SR][t] = SR(ρ)
    end
    if :κEA in observables
        data[:κEA][t] = κEA(ρ; ref=ref)
    end
    if :κ2 in observables
        data[:κ2][t] = κ2(ρ; ref=ref)
    end
    if :maxlinkdim in observables
        data[:maxlinkdim][t] = maxlinkdim(ρ)
    end
    return data
end

function circuit(L::Int, T::Int, λ::Float64, δ::Float64, q::Float64; observables=Symbol[], cutoff=1E-8, maxdim=200)
    ρ, sites = ghz(L; ref=true)

    data = Dict([s => zeros(Float64, 2T+2) for s in observables])
    data = MPS_update_data(ρ, observables, data, 1; ref=true)
    data = MPS_update_data(ρ, observables, data, 2; ref=true)

    λzz = δ*(1-λ)
    λx = δ*λ

    Xn = decoherence_layer(sites, PauliX, q, 1:L, L)
    ZZn = decoherence_layer(sites, kron(PauliZ, PauliZ), q, 1:L-1, L)
    bell_state = bell(sites)

    for t in 1:T
        ρ, _, _ = doubled_measure(ρ, kron(PauliZ,PauliZ), λzz, 1:L-1; ref=true)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        ρ = apply(ZZn, ρ)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        data = MPS_update_data(ρ, observables, data, 2t+1; ref=true)

        ρ, _, _ = doubled_measure(ρ, PauliX, λx, 1:L; ref=true)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        ρ = apply(Xn, ρ)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        data = MPS_update_data(ρ, observables, data, 2t+2; ref=true)
    end

    return ρ, data
end

function sample(L::Int, T::Int, Δ::Float64, δ::Float64, q::Float64, samples::Int; observables=Symbol[], cutoff=1E-8, maxdim=200)
    mean_data = Dict([s => zeros(Float64, 2T+2) for s in observables])
    var_data = Dict([s => zeros(Float64, 2T+2) for s in observables])
    for _ in 1:samples
        _, sample_data = circuit(L, T, Δ, δ, q; observables=observables, cutoff=cutoff, maxdim=maxdim)
        for observable in observables
            mean_data[observable] += sample_data[observable]
            var_data[observable] += sample_data[observable].^2
        end
    end
    for observable in observables
        mean_data[observable] /= samples
        var_data[observable] /= samples
    end
    return mean_data, var_data
end