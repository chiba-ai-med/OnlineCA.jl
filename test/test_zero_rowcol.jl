#####################################
# Phase 4: zero row / column handling
#
# Build a small contingency table that contains all-zero rows and all-zero
# columns, then run all three formats and assert:
#   - no NaN / Inf in any returned matrix
#   - valid_rows / valid_cols correctly identify the non-zero indices
#   - rowmass[i] == 0 for excluded rows (and zero coordinates / cos²)
#####################################

println("####### Zero rows / columns handled safely #######")

using OnlinePCA: csv2bin, mm2bin, bincoo2bin

zero_tmp = mktempdir()

# 8x6 contingency table with row 4 and column 5 set to zero.
zdata = [
    5  3  0  2  0  1
    2  4  1  0  0  3
    1  2  3  4  0  2
    0  0  0  0  0  0   # zero row
    3  1  2  0  0  4
    0  2  4  3  0  5
    4  0  1  2  0  3
    1  3  2  0  0  4
]

zN, zM = size(zdata)
@assert all(sum(zdata, dims=2)[4] .== 0)   # row 4
@assert all(sum(zdata, dims=1)[5] .== 0)   # col 5

# Write CSV / MM / BinCOO
write_csv(joinpath(zero_tmp, "Z.csv"), zdata)
mmwrite(joinpath(zero_tmp, "Z.mtx"), sparse(zdata))
zb = joinpath(zero_tmp, "Z.bincoo")
open(zb, "w") do io
    for i in 1:zN, j in 1:zM
        if zdata[i, j] != 0
            println(io, "$i $j")
        end
    end
end

csv2bin(csvfile=joinpath(zero_tmp, "Z.csv"),    binfile=joinpath(zero_tmp, "Z.zst"))
mm2bin(mmfile=joinpath(zero_tmp, "Z.mtx"),       binfile=joinpath(zero_tmp, "Z.mtx.zst"))
bincoo2bin(bincoofile=zb,                         binfile=joinpath(zero_tmp, "Z.bincoo.zst"))

zero_inputs = [
    ("ca",         joinpath(zero_tmp, "Z.zst")),
    ("sparse_ca",  joinpath(zero_tmp, "Z.mtx.zst")),
    ("bincoo_ca",  joinpath(zero_tmp, "Z.bincoo.zst")),
]

zero_funcs = Dict(
    "ca"        => ca,
    "sparse_ca" => sparse_ca,
    "bincoo_ca" => bincoo_ca,
)

# Helper: assert no NaN/Inf in every numeric matrix or vector field.
function assert_finite(out)
    for k in propertynames(out)
        v = getproperty(out, k)
        if v isa AbstractArray && eltype(v) <: AbstractFloat
            @test !any(isnan, v)
            @test !any(isinf, v)
        end
    end
end

for (name, input) in zero_inputs
    println("    -- ", name)
    f = zero_funcs[name]
    # min(N=8, M=6) is 6; dim + noversamples must be <= min(N, M).
    out = f(input=input, dim=2, noversamples=2, niter=3, chunksize=4)

    # No NaN / Inf anywhere
    assert_finite(out)

    # valid_rows / valid_cols pinpoint the surviving indices
    @test out.valid_rows == [1, 2, 3, 5, 6, 7, 8]
    @test out.valid_cols == [1, 2, 3, 4, 6]

    # Zero-mass row/col entries themselves stay at zero
    @test out.rowmass[4]  == 0
    @test out.colmass[5]  == 0
    @test all(out.rowcoord[4, :]   .== 0)
    @test all(out.colcoord[5, :]   .== 0)
    @test all(out.rowstd[4, :]     .== 0)
    @test all(out.colstd[5, :]     .== 0)
    @test all(out.rowcos2[4, :]    .== 0)
    @test all(out.colcos2[5, :]    .== 0)
end
#####################################
