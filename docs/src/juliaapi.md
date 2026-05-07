# OnlineCA.jl (Julia API)

## Correspondence Analysis (CA)
```@docs
ca(;input::AbstractString="", outdir::Union{Nothing,AbstractString}=nothing, dim::Number=3, noversamples::Number=5, niter::Number=3, chunksize::Number=0)
```

## Sparse Correspondence Analysis (SparseCA)
```@docs
sparse_ca(;input::AbstractString="", outdir::Union{Nothing,AbstractString}=nothing, dim::Number=3, noversamples::Number=5, niter::Number=3, chunksize::Number=0)
```

## BinCOO Correspondence Analysis (BinCOOCA)
```@docs
bincoo_ca(;input::AbstractString="", outdir::Union{Nothing,AbstractString}=nothing, dim::Number=3, noversamples::Number=5, niter::Number=3, chunksize::Number=0)
```
