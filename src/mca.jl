"""
Multiple Correspondence Analysis (MCA).

Given a categorical table with `q` variables, MCA is CA applied to the
N × M complete disjunctive (indicator) matrix where each row has exactly
`q` ones. We materialize the indicator as a Binary COO file on disk and
call `bincoo_ca`, so all out-of-core, sparse, randomized-SVD and
reproducibility machinery applies unchanged.

The CA spectrum on an indicator matrix is artificially inflated and
contains structurally unhelpful axes near 1/q. Phase 9 adds the standard
Benzécri / Greenacre eigenvalue corrections, selectable via the
`correction` keyword.
"""

# ---------------------------------------------------------------------------
# Indicator construction (preprocessing helpers)
# ---------------------------------------------------------------------------

"""
    categorical_to_bincoo(table, bincoofile) -> meta

Write the indicator matrix of `table` to `bincoofile` in BinCOO text
format (one `i j` line per nonzero). `table` is `N × q`; entries are
positive integer category codes. Returns a NamedTuple of metadata used
by `mca` to label categories.
"""
function categorical_to_bincoo(table::AbstractMatrix{<:Integer},
                               bincoofile::AbstractString)
    N, q = size(table)
    levels_per_var = [Int(maximum(@view(table[:, j]))) for j in 1:q]
    offsets        = vcat(0, cumsum(levels_per_var)[1:end-1])
    M              = sum(levels_per_var)

    open(bincoofile, "w") do io
        @inbounds for i in 1:N, j in 1:q
            level = table[i, j]
            level >= 1 || throw(ArgumentError("category codes must be positive integers; got $level at ($i, $j)"))
            level <= levels_per_var[j] || throw(ArgumentError("category code $level exceeds maximum for variable $j"))
            col = offsets[j] + level
            println(io, "$i $col")
        end
    end

    return (N=N, q=q, M=M,
            levels_per_var=levels_per_var,
            offsets=offsets)
end

"""
    _make_category_labels(meta, var_names) -> (categories, var_of_category)

Build the category-name vector and the category→variable index map.
"""
function _make_category_labels(meta, var_names)
    categories      = String[]
    var_of_category = Int[]
    for j in 1:meta.q
        for level in 1:meta.levels_per_var[j]
            push!(categories,      string(var_names[j], "=", level))
            push!(var_of_category, j)
        end
    end
    return categories, var_of_category
end

# ---------------------------------------------------------------------------
# Eigenvalue / inertia correction (Phase 9)
# ---------------------------------------------------------------------------

"""
    _mca_correction(sigma, q, correction) -> (adjusted_inertia, adjusted_total)

Apply the Benzécri (1979) per-eigenvalue correction. Both `:benzecri`
and `:greenacre` use the same per-axis formula

    λ*_k = ((q / (q-1)) * (λ_k - 1/q))^2  if  λ_k > 1/q,  else 0

with λ_k = σ_k². They differ in how the total adjusted inertia is
reported:

  * :benzecri  — total = Σ λ*_k                (sum of corrected axes)
  * :greenacre — total = (q/(q-1))² * Σ (λ_k - 1/q)²
                 over all axes with λ_k > 1/q (Greenacre 1993)

`correction = :none` returns the raw values unchanged.
"""
function _mca_correction(sigma::AbstractVector, q::Integer, correction::Symbol)
    raw = sigma .^ 2
    if correction === :none
        return raw, sum(raw)
    end
    correction in (:benzecri, :greenacre) ||
        throw(ArgumentError("correction must be :none, :benzecri, or :greenacre; got $correction"))

    threshold = 1 / q
    factor    = q / (q - 1)

    adjusted = similar(raw)
    @inbounds for k in eachindex(raw)
        λ = raw[k]
        adjusted[k] = λ > threshold ? (factor * (λ - threshold))^2 : zero(λ)
    end

    total = if correction === :benzecri
        sum(adjusted)
    else  # :greenacre
        s = zero(eltype(raw))
        for λ in raw
            if λ > threshold
                s += (λ - threshold)^2
            end
        end
        factor^2 * s
    end

    return adjusted, total
end

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

"""
    mca(table::AbstractMatrix{<:Integer};
        dim=3, correction=:none, var_names=nothing,
        noversamples=5, niter=3, chunksize=0,
        rng=default_rng(), seed=nothing, T=Float64,
        outdir=nothing, keep_indicator_file::Bool=false)

In-memory MCA. `table` is an `N × q` matrix where each entry is the
positive integer category code for variable `j` of observation `i`.
Internally the function

  1. materializes the indicator matrix to a temporary Bincoo + Zstandard
     file on disk,
  2. runs `bincoo_ca` on it,
  3. wraps the result with MCA metadata (`variables`, `categories`,
     `var_of_category`, `q`, `correction`) and the optionally corrected
     eigenvalues (`inertia_adjusted`, `total_inertia_adjusted`).

The temporary file is removed on return unless `keep_indicator_file` is
set to true (in which case the path is included in the result for
reuse / inspection).
"""
function mca(table::AbstractMatrix{<:Integer};
             dim::Integer=3,
             correction::Symbol=:none,
             var_names=nothing,
             noversamples::Integer=5,
             niter::Integer=3,
             chunksize::Integer=0,
             rng::AbstractRNG=default_rng(),
             seed::Union{Nothing,Integer}=nothing,
             T::Type=Float64,
             outdir::Union{Nothing,AbstractString}=nothing,
             keep_indicator_file::Bool=false)

    N, q = size(table)
    bincoofile = tempname() * ".bincoo"
    binfile    = tempname() * ".bincoo.zst"
    meta = categorical_to_bincoo(table, bincoofile)

    try
        OnlinePCA.bincoo2bin(bincoofile=bincoofile, binfile=binfile)
    finally
        rm(bincoofile; force=true)
    end

    base = bincoo_ca(input=binfile, dim=dim,
                     noversamples=noversamples, niter=niter,
                     chunksize=chunksize,
                     rng=rng, seed=seed, T=T)

    if !keep_indicator_file
        rm(binfile; force=true)
    end

    var_names_resolved = var_names === nothing ?
        [string("var", j) for j in 1:q] : collect(var_names)
    length(var_names_resolved) == q ||
        throw(ArgumentError("var_names must have length $q (got $(length(var_names_resolved)))"))

    categories, var_of_category = _make_category_labels(meta, var_names_resolved)

    inertia_adj, total_adj = _mca_correction(base.sigma, q, correction)

    extras = (variables               = var_names_resolved,
              categories              = categories,
              var_of_category         = var_of_category,
              q                       = q,
              correction              = correction,
              inertia_adjusted        = inertia_adj,
              total_inertia_adjusted  = total_adj)

    if keep_indicator_file
        extras = merge(extras, (indicator_file=binfile,))
    end

    result = merge(base, extras)

    if outdir isa String
        mkpath(outdir)
        output(outdir, result)
        # Additional MCA-specific outputs
        write_csv(joinpath(outdir, "Categories.csv"),               categories)
        write_csv(joinpath(outdir, "Variables.csv"),                var_names_resolved)
        write_csv(joinpath(outdir, "Var_of_category.csv"),          var_of_category)
        write_csv(joinpath(outdir, "Inertia_adjusted.csv"),         inertia_adj)
        write_csv(joinpath(outdir, "Total_inertia_adjusted.csv"),   total_adj)
    end

    return result
end
