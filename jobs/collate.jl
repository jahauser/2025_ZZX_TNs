# jobs/collate.jl
using JLD2
using Glob

outdir = joinpath(@__DIR__,  "../output/")
files = glob("*.jld2", outdir)

if isempty(files)
    println("No .jld2 files found in ", outdir)
    exit(0)
end

version = :PBC

# results[(L,T,lambda,delta,q)] = (sum_E, sum_E2, total_samples, obs, total_time)
# results = Dict{NTuple{5,Any}, Tuple{Dict{Symbol,Vector{ComplexF64}},
                                    # Dict{Symbol,Vector{ComplexF64}},
                                    # Int, Vector{Symbol}, Float64}}()

if version == :no_unitary
    results = Dict{NTuple{5,Any}, Tuple{Dict{Symbol,Vector{ComplexF64}},
                                    Dict{Symbol,Vector{ComplexF64}},
                                    Int, Vector{Symbol}, Float64}}()
elseif version in [:unlabeled_unitary, :full, :PBC]
    results = Dict{NTuple{7,Any}, Tuple{Dict{Symbol,Vector{ComplexF64}},
                                    Dict{Symbol,Vector{ComplexF64}},
                                    Int, Vector{Symbol}, Float64}}()
end

for f in files
    if version == :no_unitary
        @load f L T lambda delta q samples obs mean_data var_data dt
        key = (L, T, lambda, delta, q)
    elseif version == :unlabeled_unitary
        if occursin("theta", f)
            m = match(r"theta0p(\d+)_pure(true|false)", f)

            theta = parse(Float64, replace(m.captures[1], "p" => ".") |> x -> "0." * x)
            pure  = m.captures[2] == "true"
            
            @load f L T lambda delta q samples obs mean_data var_data dt
            key = (L, T, lambda, delta, q, theta, pure)
        end
    elseif version == :full
        if occursin("theta", f)
            try 
                @load f L T lambda delta q theta pure samples obs mean_data var_data dt
                key = (L, T, lambda, delta, q, theta, pure)
            catch KeyError
                @load f L T lambda delta q samples obs mean_data var_data dt
                if samples == 1
                    continue
                else
                    println(samples)
                    continue
                end
            end
        else
            println(f)
            continue
        end
    elseif version == :PBC
        if occursin("PBC", f)
            @load f L T lambda delta q theta pure samples obs mean_data var_data dt
            key = (L, T, lambda, delta, q, theta, pure)
        else
            continue
        end
    end



    # Convert per-file averages to weighted sums
    # sum_E   += samples * E[x]
    # sum_E2  += samples * E[x^2]
    if !haskey(results, key)
        sum_E  = Dict{Symbol,Vector{ComplexF64}}()
        sum_E2 = Dict{Symbol,Vector{ComplexF64}}()
        for s in obs
            sum_E[s]  = samples .* copy(mean_data[s])
            sum_E2[s] = samples .* copy(var_data[s])
        end
        results[key] = (sum_E, sum_E2, samples, obs, float(dt))
    else
        sum_E, sum_E2, acc_samples, acc_obs, acc_dt = results[key]
        @assert acc_obs == obs "Observable mismatch for $f"
        for s in obs
            sum_E[s]  .+= samples .* mean_data[s]
            sum_E2[s] .+= samples .* var_data[s]
        end
        results[key] = (sum_E, sum_E2, acc_samples + samples, acc_obs, acc_dt + dt)
    end
end

# Normalize to weighted averages: E[x] and E[x^2]
for (key, (sum_E, sum_E2, total_samples, obs, total_time)) in results
    for s in obs
        sum_E[s]  ./= total_samples
        sum_E2[s] ./= total_samples
    end
    # Overwrite with normalized values
    results[key] = (sum_E, sum_E2, total_samples, obs, total_time)
end

outfile = joinpath(outdir, "collated_results.jld2")
@save outfile results
println("Collated $(length(results)) parameter sets from $(length(files)) files -> ", outfile)
