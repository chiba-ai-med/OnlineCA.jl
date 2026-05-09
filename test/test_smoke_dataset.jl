#####################################
# External-reference correctness check on a classic CA dataset.
#
# The "smoke" dataset is the textbook example used by Greenacre (2017,
# "Correspondence Analysis in Practice", 3rd ed.) and ships with R's
# `ca` package. Eigenvalues from `ca::ca(smoke)$sv^2` are public and
# stable across implementations:
#
#   λ_1 ≈ 0.07476   (~ 87.76 % of total inertia)
#   λ_2 ≈ 0.01001   (~ 11.76 %)
#   λ_3 ≈ 0.00041   (~  0.48 %)
#   total inertia ≈ 0.08519
#
# OnlineCA's OOC implementation must reproduce these to within the
# randomized-SVD tolerance.
#####################################

println("####### Reference correctness vs Greenacre smoke dataset #######")

# 5 × 4 contingency table (Greenacre 2017 §4)
smoke = [
     4   2   3   2     # SeniorMan
     4   3   7   4     # JuniorMan
    25  10  12   4     # SeniorEmployee
    18  24  33  13     # JuniorEmployee
    10   6   7   2     # Secretary
]
@test sum(smoke) == 193   # Greenacre total

smoke_dir = mktempdir()
write_csv(joinpath(smoke_dir, "smoke.csv"), smoke)
csv2bin(csvfile=joinpath(smoke_dir, "smoke.csv"),
        binfile=joinpath(smoke_dir, "smoke.zst"))

out = ca(input=joinpath(smoke_dir, "smoke.zst"),
         dim=3, noversamples=1, niter=5, chunksize=5, seed=42)

# Published values from R's ca::ca(smoke)$sv^2
ref_eig   = [0.07476, 0.01001, 0.00041]
ref_total = 0.08519

@test isapprox(out.inertia[1],     ref_eig[1];   atol=1e-4)
@test isapprox(out.inertia[2],     ref_eig[2];   atol=1e-4)
@test isapprox(out.inertia[3],     ref_eig[3];   atol=1e-4)
@test isapprox(out.total_inertia,  ref_total;    atol=1e-4)
#####################################
