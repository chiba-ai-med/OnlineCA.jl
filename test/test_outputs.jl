#####################################
# Phase 3: rich CA outputs
#
# masses / standard coordinates / principal coordinates / contributions / cos²
# Tests are uniform across the three formats; each runs `f(input)` once with
# the same dim and validates shapes + CA-defining identities.
#####################################

println("####### Rich CA outputs (masses, std/principal coords, contributions, cos²) #######")

phase3_inputs = [
    ("ca",         joinpath(tmp, "Data.zst")),
    ("sparse_ca",  joinpath(tmp, "Data.mtx.zst")),
    ("bincoo_ca",  joinpath(tmp, "Data.bincoo.zst")),
]

phase3_funcs = Dict(
    "ca"        => ca,
    "sparse_ca" => sparse_ca,
    "bincoo_ca" => bincoo_ca,
)

const _DIM = 2

for (name, input) in phase3_inputs
    println("    -- ", name)
    f = phase3_funcs[name]
    out = f(input=input, dim=_DIM, noversamples=10, niter=3, chunksize=100)

    # ---- Shapes
    N = size(out.rowcoord, 1)
    M = size(out.colcoord, 1)

    @test size(out.rowcoord)   == (N, _DIM)
    @test size(out.colcoord)   == (M, _DIM)
    @test size(out.rowstd)     == (N, _DIM)
    @test size(out.colstd)     == (M, _DIM)
    @test size(out.rowcontrib) == (N, _DIM)
    @test size(out.colcontrib) == (M, _DIM)
    @test size(out.rowcos2)    == (N, _DIM)
    @test size(out.colcos2)    == (M, _DIM)
    @test length(out.rowmass)  == N
    @test length(out.colmass)  == M

    # ---- Mass identities
    @test all(out.rowmass .>= 0)
    @test all(out.colmass .>= 0)
    @test isapprox(sum(out.rowmass), 1.0; atol=1e-10)
    @test isapprox(sum(out.colmass), 1.0; atol=1e-10)

    # ---- Principal = standard .* sigma' (relationship)
    @test maximum(abs.(out.rowcoord .- out.rowstd .* out.sigma'))  < 1e-10
    @test maximum(abs.(out.colcoord .- out.colstd .* out.sigma'))  < 1e-10

    # ---- Contributions: each axis (column) sums to 1 because U/V are
    #      orthonormal columns. With a randomized SVD this only holds at
    #      machine precision when niter is large enough, but for the test
    #      data the sum should still be very close to 1.
    for k in 1:_DIM
        @test isapprox(sum(out.rowcontrib[:, k]), 1.0; atol=1e-3)
        @test isapprox(sum(out.colcontrib[:, k]), 1.0; atol=1e-3)
    end
    @test all(out.rowcontrib .>= 0)
    @test all(out.colcontrib .>= 0)

    # ---- cos² is bounded in [0, 1] per cell, and per-row sum over kept
    #      dims is at most 1 (becomes exactly 1 when dim covers all axes).
    @test all(out.rowcos2 .>= 0)
    @test all(out.colcos2 .>= 0)
    @test all(out.rowcos2 .<= 1 + 1e-8)
    @test all(out.colcos2 .<= 1 + 1e-8)
    @test all(sum(out.rowcos2, dims=2) .<= 1 + 1e-8)
    @test all(sum(out.colcos2, dims=2) .<= 1 + 1e-8)
end
#####################################
