function MPS_update_data(ρ::MPS, observables::Vector{Symbol}, data::Dict{Symbol,Vector{ComplexF64}}, t::Int; ref=false, max_t=0)
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

    if :ZZ in observables
        data[:ZZ][t] = boundary_ZZ_susceptibility(ρ; ref=ref)
    end
    if :ZZZZ in observables
        data[:ZZZZ][t] = boundary_ZZZZ_susceptibility(ρ; ref=ref)
    end

    if :pure_SR in observables
        data[:pure_SR][t] = pure_SR(ρ)
    end
    if :pure_κEA in observables
        data[:pure_κEA][t] = pure_κEA(ρ; ref=ref)
    end

    if max_t != 0 && max_t - t < 4
        if :terminal_order in observables
            data[:terminal_order][t] = pure_κEA(ρ; ref=ref)
        end
        if :terminal_disorder in observables
            data[:terminal_disorder][t] = disorder_EA(ρ; ref=ref)
        end
    end

    if :binder_EA in observables
        data[:binder_EA][t] = binder_EA(ρ; ref=ref)
    end

    if :four_point_EA in observables
        data[:four_point_EA][t] = four_point_EA(ρ; ref=ref)
    end
    return data
end

function circuit(L::Int, T::Int, λ::Float64, δ::Float64, q::Float64, θx::Float64, θzz::Float64; ref=true, final_perfect=false, observables=Symbol[], cutoff=1E-8, maxdim=200)
    ρ, sites = ghz(L; ref=ref)

    data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    data = MPS_update_data(ρ, observables, data, 1; ref=ref)
    data = MPS_update_data(ρ, observables, data, 2; ref=ref)

    λzz = δ*(1-λ)
    λx = δ*λ

    Xn = decoherence_layer(sites, PauliX, q, 1:L, L)
    ZZn = decoherence_layer(sites, kron(PauliZ, PauliZ), q, 1:L-1, L)

    Xc = coherent_layer(sites, PauliX, θx, 1:L, L)
    ZZc = coherent_layer(sites, kron(PauliZ, PauliZ), θzz, 1:L-1, L)
    bell_state = bell(sites)

    for t in 1:T
        ρ, _, _ = doubled_measure(ρ, kron(PauliZ,PauliZ), λzz, 1:L-1; ref=ref)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        if q > 0.0
            ρ = apply(ZZn, ρ)
            ρ /= inner(bell_state, ρ)
            truncate!(ρ; cutoff=cutoff, maxdim=maxdim)
        end

        if θzz > 0.0
            ρ = apply(ZZc, ρ)
            ρ /= inner(bell_state, ρ)
            truncate!(ρ; cutoff=cutoff, maxdim=maxdim)
        end

        data = MPS_update_data(ρ, observables, data, 2t+1; ref=ref)

        ρ, _, _ = doubled_measure(ρ, PauliX, λx, 1:L; ref=ref)
        ρ /= inner(bell_state, ρ)
        truncate!(ρ; cutoff=cutoff, maxdim=maxdim)

        if q > 0.0
            ρ = apply(Xn, ρ)
            ρ /= inner(bell_state, ρ)
            truncate!(ρ; cutoff=cutoff, maxdim=maxdim)
        end

        if θx > 0.0
            ρ = apply(Xc, ρ)
            ρ /= inner(bell_state, ρ)
            truncate!(ρ; cutoff=cutoff, maxdim=maxdim)
        end

        if t == T && final_perfect
            ρ, _, _ = doubled_measure(ρ, kron(PauliZ,PauliZ), 1.0, 1:L-1; ref=ref)
            ρ /= inner(bell_state, ρ)
            truncate!(ρ; cutoff=cutoff, maxdim=maxdim)
        end

        data = MPS_update_data(ρ, observables, data, 2t+2; ref=ref)
    end

    return ρ, data
end

function pure_circuit(L::Int, T::Int, λ::Float64, δ::Float64, θx::Float64, θzz::Float64; ref=true, PBC=false, final_perfect=false, observables=Symbol[], cutoff=1E-8, maxdim=200)
    ψ, sites = pure_ghz(L; ref=ref)

    data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    data = MPS_update_data(ψ, observables, data, 1; ref=ref, max_t=2T+2)
    data = MPS_update_data(ψ, observables, data, 2; ref=ref, max_t=2T+2)

    λzz = δ*(1-λ)
    λx = δ*λ

    Xc = singled_coherent_layer(sites, PauliX, θx, 1:L, L)

    if !PBC
        ZZc = singled_coherent_layer(sites, kron(PauliZ, PauliZ), θzz, 1:L-1, L)
    else
        ZZc = singled_coherent_layer(sites, kron(PauliZ, PauliZ), θzz, 1:L, L)
    end

    for t in 1:T
        if !PBC
            ψ, _, _ = singled_measure(ψ, kron(PauliZ,PauliZ), λzz, 1:L-1; ref=ref)
        else
            ψ, _, _ = singled_measure(ψ, kron(PauliZ,PauliZ), λzz, 1:L; ref=ref)
        end
        ψ /= norm(ψ)
        truncate!(ψ; cutoff=cutoff, maxdim=maxdim)

        if θzz > 0.0
            ψ = apply(ZZc, ψ)
            ψ /= norm(ψ)
            truncate!(ψ; cutoff=cutoff, maxdim=maxdim)
        end

        data = MPS_update_data(ψ, observables, data, 2t+1; ref=ref, max_t=2T+2)

        ψ, _, _ = singled_measure(ψ, PauliX, λx, 1:L; ref=ref)
        ψ /= norm(ψ)
        truncate!(ψ; cutoff=cutoff, maxdim=maxdim)
        # println(expect(ψ, "X"))

        if θx > 0.0
            ψ = apply(Xc, ψ)
            ψ /= norm(ψ)
            truncate!(ψ; cutoff=cutoff, maxdim=maxdim)
        end

        if t == T && final_perfect
            ψ, _, _ = singled_measure(ψ, kron(PauliZ,PauliZ), 1.0, 1:L-1; ref=ref)
            ψ /= norm(ψ)
            truncate!(ψ; cutoff=cutoff, maxdim=maxdim) 
        end

        data = MPS_update_data(ψ, observables, data, 2t+2; ref=ref, max_t=2T+2)
    end

    return ψ, data
end

function sample(L::Int, T::Int, Δ::Float64, δ::Float64, q::Float64, θx::Float64, θzz::Float64, samples::Int; ref=true, final_perfect=false, observables=Symbol[], cutoff=1E-8, maxdim=200)
    mean_data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    var_data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    for _ in 1:samples
        _, sample_data = circuit(L, T, Δ, δ, q, θx, θzz; ref=ref, final_perfect=final_perfect, observables=observables, cutoff=cutoff, maxdim=maxdim)
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

function pure_sample(L::Int, T::Int, Δ::Float64, δ::Float64, θx::Float64, θzz::Float64, samples::Int; ref=true, PBC=false, final_perfect=false, observables=Symbol[], cutoff=1E-8, maxdim=200)
    mean_data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    var_data = Dict([s => zeros(ComplexF64, 2T+2) for s in observables])
    for _ in 1:samples
        _, sample_data = pure_circuit(L, T, Δ, δ, θx, θzz; ref=ref, PBC=PBC, final_perfect=final_perfect, observables=observables, cutoff=cutoff, maxdim=maxdim)
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