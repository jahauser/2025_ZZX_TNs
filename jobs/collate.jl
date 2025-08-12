# jobs/collate.jl
using JLD2
using Glob

outdir = joinpath(@__DIR__, "..", "output")
files = glob("*.jld2", outdir)

if isempty(files)
    println("No .jld2 files found in ", outdir)
    exit(0)
end

# results[(L,T,lambda,delta,q)] = (sum_E, sum_E2, total_samples, obs, total_time)
results = Dict{NTuple{5,Any}, Tuple{Dict{Symbol,Vector{Float64}},
                                    Dict{Symbol,Vector{Float64}},
                                    Int, Vector{Symbol}, Float64}}()

for f in files
    @load f L T lambda delta q samples obs mean_data var_data dt
    key = (L, T, lambda, delta, q)

    # Convert per-file averages to weighted sums
    # sum_E   += samples * E[x]
    # sum_E2  += samples * E[x^2]
    if !haskey(results, key)
        sum_E  = Dict{Symbol,Vector{Float64}}()
        sum_E2 = Dict{Symbol,Vector{Float64}}()
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
