#####################################
println("####### CA (Dense, Julia API) #######")
tmp_ca = mktempdir()

# Use dim=2 for accuracy tests (first 2 singular values are well-separated)
out_ca = ca(input=joinpath(tmp, "Data.zst"), dim=2, noversamples=10, niter=3, chunksize=100)
F, G, sigma, inertia, total_inertia = out_ca

# Size tests
@test size(F) == (300, 2)
@test size(G) == (99, 2)
@test length(sigma) == 2
@test length(inertia) == 2
@test total_inertia isa Float64

# Singular values should be positive and decreasing
@test all(sigma .> 0)
@test issorted(sigma, rev=true)

# Compare with reference (signs may differ, so compare absolute values)
for k in 1:2
    # Singular values should be close
    @test abs(sigma[k] - ref_sigma[k]) / ref_sigma[k] < 0.05
    # Coordinates should match up to sign
    cos_sim_F = abs(dot(F[:, k], ref_F[:, k])) / (norm(F[:, k]) * norm(ref_F[:, k]))
    cos_sim_G = abs(dot(G[:, k], ref_G[:, k])) / (norm(G[:, k]) * norm(ref_G[:, k]))
    @test cos_sim_F > 0.95
    @test cos_sim_G > 0.95
end

# Total inertia should match
@test abs(total_inertia - ref_total_inertia) / ref_total_inertia < 0.05

# Inertia should be sigma^2
@test all(abs.(inertia .- sigma .^ 2) .< 1e-10)

# Test with dim=3 for size check
out_ca3 = ca(input=joinpath(tmp, "Data.zst"), dim=3, noversamples=10, niter=3, chunksize=100)
@test size(out_ca3[1]) == (300, 3)
@test size(out_ca3[2]) == (99, 3)

# Test with output directory
out_ca_dir = ca(input=joinpath(tmp, "Data.zst"), outdir=tmp_ca, dim=2, noversamples=10, niter=3, chunksize=100)
testfilesize(true,
    joinpath(tmp_ca, "Row_coordinates.csv"),
    joinpath(tmp_ca, "Col_coordinates.csv"),
    joinpath(tmp_ca, "Singular_values.csv"),
    joinpath(tmp_ca, "Inertia.csv"),
    joinpath(tmp_ca, "Total_inertia.csv"))

# Non-divisible chunk size
out_chunk_odd = ca(input=joinpath(tmp, "Data.zst"), dim=2, chunksize=51)
@test size(out_chunk_odd[1]) == (300, 2)
@test size(out_chunk_odd[2]) == (99, 2)

# Edge cases
# Single row data
single_row_data = copy(data)[1, :]'
writedlm(joinpath(tmp_ca, "SingleRow.csv"), single_row_data, ",")
csv2bin(csvfile=joinpath(tmp_ca, "SingleRow.csv"), binfile=joinpath(tmp_ca, "SingleRow.zst"))
@test_throws Exception ca(input=joinpath(tmp_ca, "SingleRow.zst"), dim=2)

# Single column data
single_col_data = copy(data)[:, 1]
writedlm(joinpath(tmp_ca, "SingleCol.csv"), single_col_data, ",")
csv2bin(csvfile=joinpath(tmp_ca, "SingleCol.csv"), binfile=joinpath(tmp_ca, "SingleCol.zst"))
@test_throws Exception ca(input=joinpath(tmp_ca, "SingleCol.zst"), dim=2)

# dim too large
@test_throws Exception ca(input=joinpath(tmp, "Data.zst"), dim=99)

# Negative dim
@test_throws Exception ca(input=joinpath(tmp, "Data.zst"), dim=-1)

# Nonexistent file
@test_throws SystemError ca(input="nonexistent_file.zst", dim=3)

#####################################

#####################################
println("####### CA (Dense, Command line) #######")
run(`$(julia) $(joinpath(bindir, "ca")) --input $(joinpath(tmp, "Data.zst")) --outdir $(tmp_ca) --dim 2 --chunksize 100`)

testfilesize(true,
    joinpath(tmp_ca, "Row_coordinates.csv"),
    joinpath(tmp_ca, "Col_coordinates.csv"),
    joinpath(tmp_ca, "Singular_values.csv"),
    joinpath(tmp_ca, "Inertia.csv"),
    joinpath(tmp_ca, "Total_inertia.csv"))
#####################################
