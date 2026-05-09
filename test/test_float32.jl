#####################################
# Phase 6: Float32 element type support
#
# Run all three formats with `T=Float32` on the same test data and assert
# that every numeric matrix in the output uses Float32, no NaN/Inf appears,
# and singular values land within ~1% of the Float64 reference.
#####################################

println("####### Float32 element type #######")

f32_inputs = [
    ("ca",         joinpath(tmp, "Data.zst")),
    ("sparse_ca",  joinpath(tmp, "Data.mtx.zst")),
    ("bincoo_ca",  joinpath(tmp, "Data.bincoo.zst")),
]

f32_funcs = Dict(
    "ca"        => ca,
    "sparse_ca" => sparse_ca,
    "bincoo_ca" => bincoo_ca,
)

for (name, input) in f32_inputs
    println("    -- ", name)
    f = f32_funcs[name]
    out_f64 = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100, seed=7, T=Float64)
    out_f32 = f(input=input, dim=2, noversamples=10, niter=3, chunksize=100, seed=7, T=Float32)

    # 1. Float32 output element types
    @test eltype(out_f32.rowcoord)   == Float32
    @test eltype(out_f32.colcoord)   == Float32
    @test eltype(out_f32.sigma)      == Float32
    @test eltype(out_f32.rowstd)     == Float32
    @test eltype(out_f32.colstd)     == Float32
    @test eltype(out_f32.rowcontrib) == Float32
    @test eltype(out_f32.colcontrib) == Float32
    @test eltype(out_f32.rowcos2)    == Float32
    @test eltype(out_f32.colcos2)    == Float32

    # 2. No NaN / Inf in any float field
    for k in propertynames(out_f32)
        v = getproperty(out_f32, k)
        if v isa AbstractArray && eltype(v) <: AbstractFloat
            @test !any(isnan, v)
            @test !any(isinf, v)
        end
    end

    # 3. Singular values close to Float64 reference (1% relative tolerance is
    #    plenty for randomized SVD with niter=3 in Float32).
    for kk in 1:length(out_f32.sigma)
        @test abs(out_f32.sigma[kk] - out_f64.sigma[kk]) / out_f64.sigma[kk] < 0.01
    end
end
#####################################
