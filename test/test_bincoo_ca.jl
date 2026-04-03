#####################################
println("####### BinCOO CA (Julia API) #######")
tmp_bincoo_ca = mktempdir()

out_bincoo_ca = bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=2, noversamples=10, niter=3, chunksize=100)
F_bc, G_bc, sigma_bc, inertia_bc, total_inertia_bc = out_bincoo_ca

# Size tests
@test size(F_bc) == (300, 2)
@test size(G_bc) == (99, 2)
@test length(sigma_bc) == 2
@test length(inertia_bc) == 2
@test total_inertia_bc isa Float64

# Singular values should be positive and decreasing
@test all(sigma_bc .> 0)
@test issorted(sigma_bc, rev=true)

# Inertia should be sigma^2
@test all(abs.(inertia_bc .- sigma_bc .^ 2) .< 1e-10)

# Total inertia should be positive
@test total_inertia_bc > 0

# Test with dim=3 for size check
out_bincoo_ca3 = bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=3, noversamples=10, niter=3, chunksize=100)
@test size(out_bincoo_ca3[1]) == (300, 3)
@test size(out_bincoo_ca3[2]) == (99, 3)

# Test with output directory
out_bincoo_ca_dir = bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), outdir=tmp_bincoo_ca, dim=2, noversamples=10, niter=3, chunksize=100)
testfilesize(true,
    joinpath(tmp_bincoo_ca, "Row_coordinates.csv"),
    joinpath(tmp_bincoo_ca, "Col_coordinates.csv"),
    joinpath(tmp_bincoo_ca, "Singular_values.csv"),
    joinpath(tmp_bincoo_ca, "Inertia.csv"),
    joinpath(tmp_bincoo_ca, "Total_inertia.csv"))

# Non-divisible chunk size
out_chunk_odd_bc = bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=2, chunksize=51)
@test size(out_chunk_odd_bc[1]) == (300, 2)
@test size(out_chunk_odd_bc[2]) == (99, 2)

# Edge cases
@test_throws Exception bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=99)
@test_throws Exception bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=-1)
@test_throws SystemError bincoo_ca(input="nonexistent_file.bincoo.zst", dim=3)

#####################################

#####################################
println("####### BinCOO CA (Command line) #######")
run(`$(julia) $(joinpath(bindir, "bincoo_ca")) --input $(joinpath(tmp, "Data.bincoo.zst")) --outdir $(tmp_bincoo_ca) --dim 2 --chunksize 100`)

testfilesize(true,
    joinpath(tmp_bincoo_ca, "Row_coordinates.csv"),
    joinpath(tmp_bincoo_ca, "Col_coordinates.csv"),
    joinpath(tmp_bincoo_ca, "Singular_values.csv"),
    joinpath(tmp_bincoo_ca, "Inertia.csv"),
    joinpath(tmp_bincoo_ca, "Total_inertia.csv"))
#####################################
