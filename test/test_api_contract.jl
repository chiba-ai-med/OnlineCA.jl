#####################################
# Phase 1: API contract snapshot tests
#
# These lock in the public-facing return-value contract of `ca`, `sparse_ca`,
# `bincoo_ca` so that refactoring (Phase 2 onward) cannot break user code that
# depends on positional destructuring, integer indexing, length, or output
# file names. Any change to these is a deliberate breaking change.
#####################################

println("####### API contract (ca / sparse_ca / bincoo_ca) #######")

contract_inputs = [
    ("ca",         joinpath(tmp, "Data.zst")),
    ("sparse_ca",  joinpath(tmp, "Data.mtx.zst")),
    ("bincoo_ca",  joinpath(tmp, "Data.bincoo.zst")),
]

contract_funcs = Dict(
    "ca"        => ca,
    "sparse_ca" => sparse_ca,
    "bincoo_ca" => bincoo_ca,
)

const EXPECTED_OUTPUT_FILES = (
    "Row_coordinates.csv",
    "Col_coordinates.csv",
    "Singular_values.csv",
    "Inertia.csv",
    "Total_inertia.csv",
)

for (name, input) in contract_inputs
    f = contract_funcs[name]
    out = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100)

    # 1. The first five fields preserve the legacy positional contract.
    #    (Phase 3 may add fields; total length grows but order of the first
    #    five is frozen.)
    @test length(out) >= 5

    # 2. Integer indexing (Tuple- or NamedTuple-compatible)
    F           = out[1]
    G           = out[2]
    sigma       = out[3]
    inertia     = out[4]
    total_iner  = out[5]
    @test F isa AbstractMatrix
    @test G isa AbstractMatrix
    @test sigma isa AbstractVector
    @test inertia isa AbstractVector
    @test total_iner isa Real

    # 3. Positional destructuring (5-LHS form keeps working even if more
    #    fields are present.)
    a, b, c, d, e = out
    @test a === F
    @test b === G
    @test c === sigma
    @test d === inertia
    @test e === total_iner

    # 4. NamedTuple field access (Phase 2 contract — additive)
    @test out isa NamedTuple
    @test out.rowcoord       === F
    @test out.colcoord       === G
    @test out.sigma          === sigma
    @test out.inertia        === inertia
    @test out.total_inertia  === total_iner
    @test propertynames(out)[1:5] == (:rowcoord, :colcoord, :sigma,
                                      :inertia, :total_inertia)

    # 5. Shapes
    @test size(F, 2) == 2
    @test size(G, 2) == 2
    @test length(sigma) == 2
    @test length(inertia) == 2

    # 6. Output file names (when outdir given)
    outdir = mktempdir()
    f(input=input, outdir=outdir, dim=2, noversamples=10, niter=3, chunksize=100)
    for fname in EXPECTED_OUTPUT_FILES
        path = joinpath(outdir, fname)
        @test isfile(path)
        @test filesize(path) != 0
    end
end
#####################################
