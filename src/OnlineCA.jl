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
    rand, AbstractRNG, MersenneTwister, default_rng
using SparseArrays:
    sparse, spzeros
using Statistics:
    mean
using OnlinePCA  # bincoo2bin used by mca's indicator preparation step

export output, parse_commandline, ca, sparse_ca, bincoo_ca,
    project_rows, project_columns,
    mca, categorical_to_bincoo

include("Utils.jl")
include("core.jl")
include("ca.jl")
include("sparse_ca.jl")
include("bincoo_ca.jl")
include("projection.jl")
include("mca.jl")

end
