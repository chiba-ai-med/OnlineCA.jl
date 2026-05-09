module OnlineCA

using ArgParse:
    ArgParseSettings, parse_args, @add_arg_table!
using CodecZstd:
    ZstdCompressorStream, ZstdDecompressorStream
using DelimitedFiles:
    writedlm, readdlm
using LinearAlgebra:
    svd, qr!, Diagonal
using ProgressMeter:
    Progress, next!
using Random:
    rand
using SparseArrays:
    sparse, spzeros
using Statistics:
    mean

export output, parse_commandline, ca, sparse_ca, bincoo_ca

include("Utils.jl")
include("core.jl")
include("ca.jl")
include("sparse_ca.jl")
include("bincoo_ca.jl")

end
