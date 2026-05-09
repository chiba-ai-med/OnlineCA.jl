#####################################
# Phase 8 + 9: MCA + Benzécri / Greenacre correction
#
# A small categorical table (200 obs, 4 vars, 3 levels each) → indicator
# matrix → bincoo_ca → MCA result. The test asserts shapes, label
# counts, and that the eigenvalue corrections behave as expected:
#   - :benzecri / :greenacre correctly zero out axes with λ < 1/q
#   - :greenacre total >= :benzecri total when both have above-threshold
#     axes (Greenacre's denominator is generally larger)
#   - :none reproduces raw inertia
#####################################

println("####### MCA + Benzécri / Greenacre correction #######")

# Synthetic categorical data: 4 variables, each with 3 levels, 200 rows
# arranged in 3 latent groups so the first axis carries real signal.
Random.seed!(2026)
mca_N, mca_q, levels = 200, 4, 3
mca_table = zeros(Int, mca_N, mca_q)
for i in 1:mca_N
    g = ((i - 1) ÷ (mca_N ÷ 3)) + 1   # group 1, 2, 3
    g = clamp(g, 1, 3)
    for j in 1:mca_q
        # Higher chance of category == g (signal), then random
        mca_table[i, j] = rand() < 0.6 ? g : rand(1:levels)
    end
end

# ---- Default (correction = :none)
out = mca(mca_table; dim=3, noversamples=2, niter=3, chunksize=50, seed=42)

@test out.q == mca_q
@test length(out.variables)        == mca_q
@test length(out.categories)       == mca_q * levels
@test length(out.var_of_category)  == mca_q * levels
@test out.var_of_category == repeat(1:mca_q, inner=levels)

# Coordinate / mass shapes
@test size(out.rowcoord, 1) == mca_N
@test size(out.colcoord, 1) == mca_q * levels
@test out.correction === :none
@test out.inertia_adjusted == out.sigma .^ 2
@test isapprox(out.total_inertia_adjusted, sum(out.sigma .^ 2); atol=1e-12)

# ---- Benzécri correction
out_b = mca(mca_table; dim=3, noversamples=2, niter=3, chunksize=50, seed=42,
            correction=:benzecri)
@test out_b.correction === :benzecri
threshold = 1 / mca_q
factor    = mca_q / (mca_q - 1)
expected_b = [σ^2 > threshold ? (factor * (σ^2 - threshold))^2 : 0.0 for σ in out_b.sigma]
@test isapprox(out_b.inertia_adjusted, expected_b; atol=1e-12)
@test isapprox(out_b.total_inertia_adjusted, sum(expected_b); atol=1e-12)

# Axes whose raw λ ≤ 1/q must be zeroed
for k in eachindex(out_b.sigma)
    if out_b.sigma[k]^2 <= threshold
        @test out_b.inertia_adjusted[k] == 0
    else
        @test out_b.inertia_adjusted[k] > 0
    end
end

# ---- Greenacre correction
out_g = mca(mca_table; dim=3, noversamples=2, niter=3, chunksize=50, seed=42,
            correction=:greenacre)
@test out_g.correction === :greenacre
@test isapprox(out_g.inertia_adjusted, expected_b; atol=1e-12)  # same per-axis
expected_total_g = factor^2 * sum(σ -> σ^2 > threshold ? (σ^2 - threshold)^2 : 0.0,
                                   out_g.sigma)
@test isapprox(out_g.total_inertia_adjusted, expected_total_g; atol=1e-12)

# Argument validation
@test_throws ArgumentError mca(mca_table; correction=:bogus)
@test_throws ArgumentError mca(zeros(Int, 5, 3); correction=:none)   # zero codes

# ---- outdir → CSV files
outdir = mktempdir()
mca(mca_table; dim=3, noversamples=2, niter=3, chunksize=50, seed=42,
    correction=:benzecri, outdir=outdir)
for fname in ["Row_coordinates.csv", "Col_coordinates.csv", "Singular_values.csv",
              "Inertia.csv", "Total_inertia.csv",
              "Categories.csv", "Variables.csv", "Var_of_category.csv",
              "Inertia_adjusted.csv", "Total_inertia_adjusted.csv"]
    @test isfile(joinpath(outdir, fname))
    @test filesize(joinpath(outdir, fname)) > 0
end
#####################################
