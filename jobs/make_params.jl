const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
format_line(L, T, lambda, delta, q, samples) =
    "--L $L --T $T --lambda $lambda --delta $delta --q $q --samples $samples"

"Write or append parameter lines to jobs/params.txt."
function write_params(lines::Vector{String}; append::Bool=false)
    open(PARAMS_FILE, append ? "a" : "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    action = append ? "Appended" : "Wrote"
    println("$action $(length(lines)) lines to $(PARAMS_FILE)")
end

# -------------------------------
# Define your sweeps here
# -------------------------------
Ls       = union([2], 8:8:40)
lambdas  = 0.0:0.1:1.0
deltas   = [0.7]
qs       = [0.1]
samples  = Dict{Int,Int}(
    2 => 1,
    8 => 1000,
    16 => 100,
    24 => 100,
    32 => 50,
    40 => 40,
)
repeats  = Dict{Int,Int}(
    2 => 1,
    8 => 1,
    16 => 10,
    24 => 10,
    32 => 20,
    40 => 25,
)
append = false  # change to false to overwrite

lines = String[]
for L in Ls, λ in lambdas, δ in deltas, q in qs
    for _ in 1:repeats[L]
        push!(lines, format_line(L, 2L+2, λ, δ, q, samples[L]))
    end
end

# Change append=true to add to existing file
write_params(lines; append=append)