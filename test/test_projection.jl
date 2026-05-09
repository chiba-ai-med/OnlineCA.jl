#####################################
# Phase 7: supplementary projection
#
# The supplementary projection formula reduces to the active formula when the
# "supplementary" row IS one of the active rows. So projecting active rows /
# columns back through `project_rows` / `project_columns` must reproduce the
# corresponding principal coordinates within randomized-SVD tolerance.
#####################################

println("####### Supplementary projection (row/col round-trip) #######")

# Use the existing 300×99 test data; pick a few active rows/cols and project.
ref3 = reference_ca(data, 3)

# ---- Row projection round-trip
proj_F = project_rows(ref3, Float64.(data[1:5, :]))
@test size(proj_F) == (5, 3)
for i in 1:5, k in 1:3
    @test isapprox(proj_F[i, k], ref3.rowcoord[i, k]; atol=1e-9)
end

# ---- Column projection round-trip
proj_G = project_columns(ref3, Float64.(data[:, 1:5]))
@test size(proj_G) == (5, 3)
for j in 1:5, k in 1:3
    @test isapprox(proj_G[j, k], ref3.colcoord[j, k]; atol=1e-9)
end

# ---- Result from the OOC ca() call: the projection function must accept
#      its NamedTuple and return the right shape / no NaN. We do not
#      assert numerical equality because U, V, σ from randomized SVD do
#      not form an exact factorization and self-consistency between
#      `rowcoord` and `colstd` only holds approximately.
out_oc = ca(input=joinpath(tmp, "Data.zst"), dim=3, noversamples=10,
            niter=3, chunksize=100, seed=42)
proj_F_oc = project_rows(out_oc,    Float64.(data[1:5, :]))
proj_G_oc = project_columns(out_oc, Float64.(data[:, 1:5]))
@test size(proj_F_oc) == (5, 3)
@test size(proj_G_oc) == (5, 3)
@test !any(isnan, proj_F_oc) && !any(isinf, proj_F_oc)
@test !any(isnan, proj_G_oc) && !any(isinf, proj_G_oc)

# ---- Argument validation
@test_throws DimensionMismatch project_rows(ref3, ones(2, 50))   # wrong M
@test_throws DimensionMismatch project_columns(ref3, ones(50, 2))  # wrong N
@test_throws ArgumentError    project_rows(ref3, zeros(2, 99))    # zero row totals
@test_throws ArgumentError    project_columns(ref3, zeros(300, 2))
#####################################
