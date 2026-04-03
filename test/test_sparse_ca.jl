#####################################
println("####### Sparse CA (Julia API) #######")
tmp_sparse_ca = mktempdir()

out_sparse_ca = sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=2, noversamples=10, niter=3, chunksize=100)
F_sp, G_sp, sigma_sp, inertia_sp, total_inertia_sp = out_sparse_ca

# Size tests
@test size(F_sp) == (300, 2)
@test size(G_sp) == (99, 2)
@test length(sigma_sp) == 2
@test length(inertia_sp) == 2
@test total_inertia_sp isa Float64

# Singular values should be positive and decreasing
@test all(sigma_sp .> 0)
@test issorted(sigma_sp, rev=true)

# Compare with reference
for k in 1:2
    @test abs(sigma_sp[k] - ref_sigma[k]) / ref_sigma[k] < 0.05
    cos_sim_F = abs(dot(F_sp[:, k], ref_F[:, k])) / (norm(F_sp[:, k]) * norm(ref_F[:, k]))
    cos_sim_G = abs(dot(G_sp[:, k], ref_G[:, k])) / (norm(G_sp[:, k]) * norm(ref_G[:, k]))
    @test cos_sim_F > 0.95
    @test cos_sim_G > 0.95
end

# Total inertia should match
@test abs(total_inertia_sp - ref_total_inertia) / ref_total_inertia < 0.05

# Test with dim=3 for size check
out_sparse_ca3 = sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=3, noversamples=10, niter=3, chunksize=100)
@test size(out_sparse_ca3[1]) == (300, 3)
@test size(out_sparse_ca3[2]) == (99, 3)

# Test with output directory
out_sparse_ca_dir = sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), outdir=tmp_sparse_ca, dim=2, noversamples=10, niter=3, chunksize=100)
testfilesize(true,
    joinpath(tmp_sparse_ca, "Row_coordinates.csv"),
    joinpath(tmp_sparse_ca, "Col_coordinates.csv"),
    joinpath(tmp_sparse_ca, "Singular_values.csv"),
    joinpath(tmp_sparse_ca, "Inertia.csv"),
    joinpath(tmp_sparse_ca, "Total_inertia.csv"))

# Non-divisible chunk size
out_chunk_odd_sp = sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=2, chunksize=51)
@test size(out_chunk_odd_sp[1]) == (300, 2)
@test size(out_chunk_odd_sp[2]) == (99, 2)

# Edge cases
@test_throws Exception sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=99)
@test_throws Exception sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=-1)
@test_throws SystemError sparse_ca(input="nonexistent_file.mtx.zst", dim=3)

#####################################

#####################################
println("####### Sparse CA (Command line) #######")
run(`$(julia) $(joinpath(bindir, "sparse_ca")) --input $(joinpath(tmp, "Data.mtx.zst")) --outdir $(tmp_sparse_ca) --dim 2 --chunksize 100`)

testfilesize(true,
    joinpath(tmp_sparse_ca, "Row_coordinates.csv"),
    joinpath(tmp_sparse_ca, "Col_coordinates.csv"),
    joinpath(tmp_sparse_ca, "Singular_values.csv"),
    joinpath(tmp_sparse_ca, "Inertia.csv"),
    joinpath(tmp_sparse_ca, "Total_inertia.csv"))
#####################################
