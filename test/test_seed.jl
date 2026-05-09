#####################################
# Phase 5: seed / rng reproducibility
#
# Same `seed` must produce bitwise-identical numeric outputs across runs.
# Different `seed` should produce numerically distinct outputs (the
# subspace spans should overlap, so cosine similarity is high, but the
# raw matrices should differ).
#####################################

println("####### Reproducibility (seed / rng) #######")

seed_inputs = [
    ("ca",         joinpath(tmp, "Data.zst")),
    ("sparse_ca",  joinpath(tmp, "Data.mtx.zst")),
    ("bincoo_ca",  joinpath(tmp, "Data.bincoo.zst")),
]

seed_funcs = Dict(
    "ca"        => ca,
    "sparse_ca" => sparse_ca,
    "bincoo_ca" => bincoo_ca,
)

for (name, input) in seed_inputs
    println("    -- ", name)
    f = seed_funcs[name]

    # Same seed → identical results bit-for-bit
    a = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100, seed=42)
    b = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100, seed=42)
    @test a.rowcoord == b.rowcoord
    @test a.colcoord == b.colcoord
    @test a.sigma    == b.sigma
    @test a.rowstd   == b.rowstd
    @test a.colstd   == b.colstd
    @test a.rowcontrib == b.rowcontrib

    # Different seed → at least some entries differ
    c = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100, seed=99)
    @test a.rowcoord != c.rowcoord

    # Passing a fresh MersenneTwister directly works the same way
    using Random
    d = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100,
          rng=MersenneTwister(42))
    @test d.rowcoord == a.rowcoord
end
#####################################
