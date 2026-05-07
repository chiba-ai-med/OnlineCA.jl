# OnlineCA.jl
Online Correspondence Analysis

[![Build Status](https://github.com/chiba-ai-med/OnlineCA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/chiba-ai-med/OnlineCA.jl/actions/workflows/CI.yml?query=branch%3Amain)

## 📚 Documentation
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://chiba-ai-med.github.io/OnlineCA.jl/dev)

## Description
OnlineCA.jl performs out-of-core Correspondence Analysis (CA) for extremely large scale matrix without loading the whole data on the memory space.

CA is a statistical method to analyze a contingency table by decomposing the standardized residual matrix

```
S = D_r^{-1/2} (P - r_p c_p') D_c^{-1/2}
```

where `P = X / n` is the relative-frequency matrix, `r_p = P 1` is the row-mass vector, `c_p = P' 1` is the column-mass vector, and `D_r = diag(r_p)`, `D_c = diag(c_p)` are the corresponding diagonal mass matrices. The matrix `S` is never explicitly formed. Instead, matrix–vector products `S * v` and `S' * u` are computed by streaming through the data, and the top singular triplets are extracted with a Halko-style randomized SVD.

__Note: The input matrix is supposed to be a non-negative count matrix (contingency table).__

## Algorithms
- Halko-style Randomized SVD on the implicit standardized residual matrix : [Halko, N. et al., 2011](https://arxiv.org/abs/0909.4061), [Halko, N. et al., 2011](https://epubs.siam.org/doi/abs/10.1137/100804139)

## Installation
```julia
# push the key "]" and type the following command.
(@julia) pkg> add https://github.com/rikenbit/OnlinePCA.jl
(@julia) pkg> add https://github.com/chiba-ai-med/OnlineCA.jl
(@julia) pkg> add PlotlyJS
# After that, push Ctrl + C to leave from Pkg REPL mode
```

## Basic API usage
### Preprocess of CSV
```julia
using OnlinePCA
using OnlinePCA: write_csv
using OnlineCA
using Distributions
using DelimitedFiles
using SparseArrays
using MatrixMarket

# CSV (Input data is supposed to be non-negative integer count)
tmp = mktempdir()
data = Int64.(ceil.(rand(NegativeBinomial(1, 0.5), 300, 99)))
data[1:50, 1:33] .= 100*data[1:50, 1:33]
data[51:100, 34:66] .= 100*data[51:100, 34:66]
data[101:150, 67:99] .= 100*data[101:150, 67:99]
write_csv(joinpath(tmp, "Data.csv"), data)

# Matrix Market (MM)
mmwrite(joinpath(tmp, "Data.mtx"), sparse(data))

# Binary COO (BinCOO)
data2 = zeros(Int, 300, 99)
data2[1:50, 1:33] .= 1
data2[51:100, 34:66] .= 1
data2[101:150, 67:99] .= 1
data2[151:300, :] .= 1

bincoofile = joinpath(tmp, "Data.bincoo")
open(bincoofile, "w") do io
    for i in 1:size(data2, 1)
        for j in 1:size(data2, 2)
            if data2[i, j] != 0
                println(io, "$i $j")
            end
        end
    end
end

# Binarization (Zstandard)
csv2bin(csvfile=joinpath(tmp, "Data.csv"), binfile=joinpath(tmp, "Data.zst"))

# Sparsification (Zstandard + MM format)
mm2bin(mmfile=joinpath(tmp, "Data.mtx"), binfile=joinpath(tmp, "Data.mtx.zst"))

# Binarization (BinCOO + Zstandard)
bincoo2bin(bincoofile=bincoofile, binfile=joinpath(tmp, "Data.bincoo.zst"))
```

### Setting for plot
```julia
using DataFrames
using PlotlyJS

function subplots(out_ca, group)
    # data frame
    data_left = DataFrame(ca1=out_ca[1][:,1], ca2=out_ca[1][:,2], group=group)
    data_right = DataFrame(ca2=out_ca[1][:,2], ca3=out_ca[1][:,3], group=group)
    # plot
    p_left = Plot(data_left, x=:ca1, y=:ca2, mode="markers", marker_size=10, group=:group)
    p_right = Plot(data_right, x=:ca2, y=:ca3, mode="markers", marker_size=10,
    group=:group, showlegend=false)
    p_left.data[1]["marker_color"] = "red"
    p_left.data[2]["marker_color"] = "blue"
    p_left.data[3]["marker_color"] = "green"
    p_right.data[1]["marker_color"] = "red"
    p_right.data[2]["marker_color"] = "blue"
    p_right.data[3]["marker_color"] = "green"
    p_left.data[1]["name"] = "group1"
    p_left.data[2]["name"] = "group2"
    p_left.data[3]["name"] = "group3"
    p_left.layout["title"] = "Component 1 vs Component 2"
    p_right.layout["title"] = "Component 2 vs Component 3"
    p_left.layout["xaxis_title"] = "ca-1"
    p_left.layout["yaxis_title"] = "ca-2"
    p_right.layout["xaxis_title"] = "ca-2"
    p_right.layout["yaxis_title"] = "ca-3"
    plot([p_left p_right])
end

group=vcat(repeat(["group1"],inner=100), repeat(["group2"],inner=100), repeat(["group3"],inner=100))
```

### Correspondence Analysis (CA)
```julia
out_ca = ca(input=joinpath(tmp, "Data.zst"), dim=3, noversamples=5, niter=3, chunksize=100)

subplots(out_ca, group)
```
![CA](./docs/src/figure/ca.png)

### Sparse Correspondence Analysis (Sparse-CA)
```julia
out_sparse_ca = sparse_ca(input=joinpath(tmp, "Data.mtx.zst"), dim=3, noversamples=5, niter=3, chunksize=100)

subplots(out_sparse_ca, group)
```
![SPARSE_CA](./docs/src/figure/sparse_ca.png)

### BinCOO Correspondence Analysis (BinCOO-CA)
```julia
out_bincoo_ca = bincoo_ca(input=joinpath(tmp, "Data.bincoo.zst"), dim=3, noversamples=5, niter=3, chunksize=100)

subplots(out_bincoo_ca, group)
```
![BinCOO_CA](./docs/src/figure/bincoo_ca.png)

## Output files
Each CA function returns a tuple `(F, G, σ, Inertia, TotalInertia)` and, when `outdir` is set, writes the following CSV files to `outdir`:

- `Row_coordinates.csv` : Row principal coordinates `F` (N × dim)
- `Col_coordinates.csv` : Column principal coordinates `G` (M × dim)
- `Singular_values.csv` : Singular values `σ` (dim,)
- `Inertia.csv` : Inertia explained by each dimension (`σ.^2`)
- `Total_inertia.csv` : Total inertia (chi-squared / n)

## Command line usage
The type of input file is assumed to be CSV, MM, or BinCOO format, and is preprocessed to a Zstandard-compressed binary file by `csv2bin`, `mm2bin`, or `bincoo2bin` in the `OnlinePCA` package. The binary file is specified as the input of the CA functions in `OnlineCA`. All CA functions can be performed as command line tools with the same parameter names like below.

```bash
# CSV → Julia Binary
julia YOUR_HOME_DIR/.julia/v0.x/OnlinePCA/bin/csv2bin \
    --csvfile Data.csv --binfile Data.zst

# MM → Julia Binary
julia YOUR_HOME_DIR/.julia/v0.x/OnlinePCA/bin/mm2bin \
    --mmfile Data.mtx --binfile Data.mtx.zst

# BinCOO → Julia Binary
julia YOUR_HOME_DIR/.julia/v0.x/OnlinePCA/bin/bincoo2bin \
    --bincoofile Data.bincoo --binfile Data.bincoo.zst

# CA (Dense)
julia YOUR_HOME_DIR/.julia/v0.x/OnlineCA/bin/ca \
    --input Data.zst --dim 3 \
    --noversamples 5 --niter 3 --chunksize 100

# Sparse-CA
julia YOUR_HOME_DIR/.julia/v0.x/OnlineCA/bin/sparse_ca \
    --input Data.mtx.zst --dim 3 \
    --noversamples 5 --niter 3 --chunksize 100

# BinCOO-CA
julia YOUR_HOME_DIR/.julia/v0.x/OnlineCA/bin/bincoo_ca \
    --input Data.bincoo.zst --dim 3 \
    --noversamples 5 --niter 3 --chunksize 100
```

## Contributing

If you have suggestions for how `OnlineCA.jl` could be improved, or want to report a bug, open an issue! We'd love all and any contributions.

## Author
- Koki Tsuyuzaki
