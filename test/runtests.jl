using OnlineCA
using OnlinePCA
using OnlinePCA: read_csv, write_csv
using Test
using DelimitedFiles
using Statistics
using Distributions
using SparseArrays
using MatrixMarket
using Random
using LinearAlgebra

# Setting
tmp = mktempdir()
julia = joinpath(Sys.BINDIR, "julia")
bindir = joinpath(dirname(pathof(OnlineCA)), "..", "bin")

function testfilesize(remove::Bool, x...)
    for n = 1:length(x)
        @test filesize(x[n]) != 0
        if remove
            rm(x[n])
        end
    end
end

# Generate test data: a contingency table with block structure
Random.seed!(123456)
data = Int64.(ceil.(rand(NegativeBinomial(1, 0.5), 300, 99)))
data[1:100, 1:33] .= 100 * data[1:100, 1:33]
data[101:200, 34:66] .= 100 * data[101:200, 34:66]
data[201:300, 67:99] .= 100 * data[201:300, 67:99]

# CSV
write_csv(joinpath(tmp, "Data.csv"), data)

# Matrix Market (MM)
sparse_data = sparse(data)
mmwrite(joinpath(tmp, "Data.mtx"), sparse_data)

# Binary COO (BinCOO)
bincoofile = joinpath(tmp, "Data.bincoo")
open(bincoofile, "w") do io
    for i in 1:size(data, 1)
        for j in 1:size(data, 2)
            if data[i, j] != 0
                println(io, "$i $j")
            end
        end
    end
end

# CSV => Zstandard
csv2bin(csvfile=joinpath(tmp, "Data.csv"),
    binfile=joinpath(tmp, "Data.zst"))

# MM => Zstandard
mm2bin(mmfile=joinpath(tmp, "Data.mtx"),
    binfile=joinpath(tmp, "Data.mtx.zst"))

# Binarization (BinCOO + Zstandard)
# Use OnlinePCA's bincoo2bin if available, otherwise OnlineNMF
bincoo2bin(bincoofile=bincoofile, binfile=joinpath(tmp, "Data.bincoo.zst"))

# Reference CA: compute exact CA in-memory for comparison.
# Returns the same NamedTuple shape as the OOC implementation (Phase 3+).
function reference_ca(X, dim)
    X = Float64.(X)
    n_total = sum(X)
    P = X ./ n_total
    r = vec(sum(P, dims=2))   # row masses
    c = vec(sum(P, dims=1))   # column masses
    Dr_inv_sqrt = Diagonal(1.0 ./ sqrt.(r))
    Dc_inv_sqrt = Diagonal(1.0 ./ sqrt.(c))
    S = Dr_inv_sqrt * (P - r * c') * Dc_inv_sqrt
    U_s, sigma_s, V_s = svd(S)

    sigma_dim = sigma_s[1:dim]
    U_dim     = U_s[:, 1:dim]
    V_dim     = V_s[:, 1:dim]

    rowstd     = Dr_inv_sqrt * U_dim
    colstd     = Dc_inv_sqrt * V_dim
    rowcoord   = rowstd .* sigma_dim'
    colcoord   = colstd .* sigma_dim'
    rowcontrib = U_dim .^ 2
    colcontrib = V_dim .^ 2

    rowdist2 = vec(sum(S .^ 2, dims=2))   # Σ_j S_ij²
    coldist2 = vec(sum(S .^ 2, dims=1))   # Σ_i S_ij²
    rowcos2  = (rowcoord .^ 2) ./ rowdist2
    colcos2  = (colcoord .^ 2) ./ coldist2

    return (rowcoord       = rowcoord,
            colcoord       = colcoord,
            sigma          = sigma_dim,
            inertia        = sigma_dim .^ 2,
            total_inertia  = sum(S .^ 2),
            rowstd         = rowstd,
            colstd         = colstd,
            rowmass        = r,
            colmass        = c,
            rowcontrib     = rowcontrib,
            colcontrib     = colcontrib,
            rowcos2        = rowcos2,
            colcos2        = colcos2)
end

ref = reference_ca(data, 3)
ref_F, ref_G, ref_sigma, ref_inertia, ref_total_inertia =
    ref.rowcoord, ref.colcoord, ref.sigma, ref.inertia, ref.total_inertia

# Tests
println("Running all tests...")

include("test_api_contract.jl")
include("test_outputs.jl")
include("test_zero_rowcol.jl")
include("test_seed.jl")
include("test_float32.jl")
include("test_ca.jl")
include("test_sparse_ca.jl")
include("test_bincoo_ca.jl")

println("All tests completed.")
