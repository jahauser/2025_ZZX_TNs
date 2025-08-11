using JLD2
using Glob

# Load all .jld2 files in output/
outdir = joinpath(@__DIR__, "..", "output")
files = glob("*.jld2", outdir)

# Collated results
# Dict: (L,T,lambda,delta,q) => (mean_data, var_data, total_samples, obs)
results = Dict{NTuple{5,Any}, Any}()

for f in files
    @load f L T lambda delta q samples obs mean_data var_data dt
    key = (L, T, lambda, delta, q)

    if !haskey(results, key)
        # Deepcopy so we don't mutate original arrays from the first file
        results[key] = (deepcopy(mean_data), deepcopy(var_data), samples, obs, dt)
    else
        acc_mean, acc_var, acc_samples, acc_obs, acc_dt = results[key]
        @assert acc_obs == obs "Observable mismatch for $f"
        for s in obs
            acc_mean[s] .+= mean_data[s] .* samples
            acc_var[s]  .+= var_data[s]  .* samples
        end
        results[key] = (acc_mean, acc_var, acc_samples + samples, acc_obs, acc_dt + dt)
    end
end

# Normalize weighted sums to averages
for (key, (m, v, total_samples, obs, total_time)) in results
    for s in obs
        m[s] ./= total_samples
        v[s] ./= total_samples
    end
    results[key] = (m, v, total_samples, obs, total_time)
end

outfile = joinpath(outdir, "collated_results.jld2")
@save outfile results
