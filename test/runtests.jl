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

# Reference CA: compute exact CA in-memory for comparison
function reference_ca(X, dim)
    X = Float64.(X)
    n_total = sum(X)
    P = X ./ n_total
    r = vec(sum(P, dims=2))  # row masses
    c = vec(sum(P, dims=1))  # column masses
    # Standardized residual matrix
    Dr_inv_sqrt = Diagonal(1.0 ./ sqrt.(r))
    Dc_inv_sqrt = Diagonal(1.0 ./ sqrt.(c))
    S = Dr_inv_sqrt * (P - r * c') * Dc_inv_sqrt
    # SVD
    U_s, sigma_s, V_s = svd(S)
    # Row and column coordinates (principal coordinates)
    F = zeros(size(X, 1), dim)
    G = zeros(size(X, 2), dim)
    for k in 1:dim
        F[:, k] = sigma_s[k] .* U_s[:, k]
        G[:, k] = sigma_s[k] .* V_s[:, k]
    end
    # Total inertia
    total_inertia = sum(S .^ 2)
    return F, G, sigma_s[1:dim], sigma_s[1:dim] .^ 2, total_inertia
end

ref_F, ref_G, ref_sigma, ref_inertia, ref_total_inertia = reference_ca(data, 3)

# Tests
println("Running all tests...")

include("test_ca.jl")
include("test_sparse_ca.jl")
include("test_bincoo_ca.jl")

println("All tests completed.")
