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
include(joinpath(ROOT, "src", "singled_measurements.jl"))
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
        "--theta"
            help = "Parameter theta"
            arg_type = Float64
            required = true
        "--samples"
            help = "Number of samples (averaging inside sample())"
            arg_type = Int
            default = 1
        "--pure"
            arg_type = Bool
            default = false
        "--PBC"
            arg_type = Bool
            default = false
    end
    return s
end

# Filename-safe run tag
function tag(; L, T, lambda, delta, q, theta, pure)
    f(x) = replace(@sprintf("%.3f", x), "." => "p")
    return "L$(L)_T$(T)_lambda$(f(lambda))_delta$(f(delta))_q$(f(q))_theta$(f(theta))_pure$pure"
end

function main(args)
    opts   = parse_args(args, build_parser())
    L      = opts["L"]
    T      = opts["T"]
    lambda = opts["lambda"]
    delta  = opts["delta"]
    q      = opts["q"]
    theta = opts["theta"]
    samples = opts["samples"]
    pure = opts["pure"]
    PBC = opts["PBC"]

    # Hardcoded observables
    if pure
        if PBC
            obs = [:terminal_order, :terminal_disorder, :pure_SR]
        else
            obs = [:pure_SR, :pure_κEA]
        end
    else
        obs = [:Ic, :SR, :κEA, :κ2, :maxlinkdim]
    end

    # Hardcoded output dir for easy local testing
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)
    
    _, _ = sample(2, 2, 0.1, 0.1, 0.1, 0.1, 0.1, 1; observables=obs, cutoff=1e-8, maxdim=200)

    t1 = time()
    if pure
        mean_data, var_data = pure_sample(L, T, lambda, delta, theta, theta, samples; PBC=PBC, ref=true, final_perfect=true, observables=obs, cutoff=1e-8, maxdim=200)
    else
        mean_data, var_data = sample(L, T, lambda, delta, q, theta, theta, samples; ref=true, final_perfect=true, observables=obs, cutoff=1e-8, maxdim=200)
    end
    dt = time()-t1

    tagstr = tag(L=L, T=T, lambda=lambda, delta=delta, q=q, theta=theta, pure=pure)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))

    if !PBC
        fname  = joinpath(outdir, "sample_$(tagstr)_$(timestr)_$(randtag).jld2")
    else
        fname  = joinpath(outdir, "sample_$(tagstr)_PBC_$(timestr)_$(randtag).jld2")
    end

    @save fname L T lambda delta q theta pure samples obs mean_data var_data dt
end

main(ARGS)
