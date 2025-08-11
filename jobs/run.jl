using ArgParse
using Printf
using Dates
using JLD2
using ITensors
using ITensorMPS
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "main.jl"))
include(joinpath(ROOT, "src", "observables.jl"))
include(joinpath(ROOT, "src", "states.jl"))
include(joinpath(ROOT, "src", "doubled_measurements.jl"))
include(joinpath(ROOT, "src", "circuits.jl"))

function build_parser()
    s = ArgParseSettings(; description = "Run sample() and save results as JLD2.")
    @add_arg_table s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--T"
            help = "Time steps"
            arg_type = Int
            required = true
        "--lambda"
            help = "Parameter lambda"
            arg_type = Float64
            required = true
        "--delta"
            help = "Parameter delta"
            arg_type = Float64
            required = true
        "--q"
            help = "Parameter q"
            arg_type = Float64
            required = true
        "--samples"
            help = "Number of samples (averaging inside sample())"
            arg_type = Int
            default = 1
    end
    return s
end

# Filename-safe run tag
function tag(; L, T, lambda, delta, q)
    f(x) = replace(@sprintf("%.3f", x), "." => "p")
    return "L$(L)_T$(T)_lambda$(f(lambda))_delta$(f(delta))_q$(f(q))"
end

function main(args)
    opts   = parse_args(args, build_parser())
    L      = opts["L"]
    T      = opts["T"]
    lambda = opts["lambda"]
    delta  = opts["delta"]
    q      = opts["q"]
    samples = opts["samples"]

    # Hardcoded observables
    obs = [:Ic, :SR, :κEA, :κ2, :maxlinkdim]

    # Hardcoded output dir for easy local testing
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)
    
    _, _ = sample(2, 2, 0.1, 0.1, 0.1, 1; observables=obs, cutoff=1e-8, maxdim=200)

    t1 = time()
    mean_data, var_data = sample(L, T, lambda, delta, q, samples; observables=obs, cutoff=1e-8, maxdim=200)
    dt = time()-t1

    tagstr = tag(L=L, T=T, lambda=lambda, delta=delta, q=q)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))
    fname  = joinpath(outdir, "sample_$(tagstr)_$(timestr)_$(randtag).jld2")

    @save fname L T lambda delta q samples obs mean_data var_data dt
end

main(ARGS)
