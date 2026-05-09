"""
Common, format-agnostic helpers shared by `ca`, `sparse_ca`, `bincoo_ca`.

The math (Greenacre 2017, "Correspondence Analysis in Practice"):

S = D_r^{-1/2} (P - r_p c_p^T) D_c^{-1/2} = U Σ V^T  (SVD)

  rowmass    : r_p              = rowsums  / total           (length N)
  colmass    : c_p              = colsums  / total           (length M)
  rowstd     : D_r^{-1/2} U                                  (N × dim)
  colstd     : D_c^{-1/2} V                                  (M × dim)
  rowcoord   : D_r^{-1/2} U Σ   = rowstd .* σ'               (N × dim)
  colcoord   : D_c^{-1/2} V Σ   = colstd .* σ'               (M × dim)
  rowcontrib : U .^ 2           ; columns sum to 1            (N × dim)
  colcontrib : V .^ 2           ; columns sum to 1            (M × dim)
  rowcos2    : (σ_k u_ik)^2 / Σ_j S_ij^2                     (N × dim)
  colcos2    : (σ_k v_jk)^2 / Σ_i S_ij^2                     (M × dim)
"""

"""
Resolve an effective RNG. When `seed === nothing`, return `rng` unchanged so
callers can thread their own state. When `seed` is given, build a fresh
`MersenneTwister(seed)` so the run is reproducible regardless of any external
RNG state.
"""
_make_rng(rng::AbstractRNG, ::Nothing) = rng
_make_rng(::AbstractRNG, seed::Integer) = MersenneTwister(seed)

# Defensive 1/sqrt: returns 0 for non-positive masses (Phase 4 will filter).
@inline _safe_inv_sqrt(x::Real) = x > 0 ? 1 / sqrt(x) : zero(x)

"Per-row inverse-sqrt mass in element type `T` (length N)."
function _inv_sqrt_mass(::Type{T}, sums::AbstractVector, total::Real) where T
    out = Vector{T}(undef, length(sums))
    @inbounds for i in eachindex(sums)
        m = sums[i] / total
        out[i] = m > 0 ? T(1 / sqrt(m)) : zero(T)
    end
    return out
end

"Standard coordinates: D^{-1/2} * U_or_V (size N×dim or M×dim). Preserves eltype(U)."
function _standard_coords(U::AbstractMatrix, inv_sqrt_mass::AbstractVector)
    N, k = size(U)
    @assert length(inv_sqrt_mass) == N
    out = similar(U)
    @inbounds for j in 1:k, i in 1:N
        out[i, j] = inv_sqrt_mass[i] * U[i, j]
    end
    return out
end

"Principal coordinates: standard .* sigma' (broadcasting σ along columns). Preserves eltype."
function _principal_coords(std::AbstractMatrix, sigma::AbstractVector)
    N, k = size(std)
    @assert length(sigma) == k
    out = similar(std)
    @inbounds for j in 1:k, i in 1:N
        out[i, j] = std[i, j] * sigma[j]
    end
    return out
end

"Row/column contributions: U.^2 (each column sums to 1 because U is orthonormal)."
function _contributions(U::AbstractMatrix)
    return U .^ 2
end

"""
Squared cosines of the angle between row/col profile and each axis.
   cos²_{ik} = (σ_k * U_{ik})^2 / dist2_i

`dist2` is the per-row (or per-col) sum of S^2 entries. Zero masses propagate
to cos²=0. Preserves eltype(U).
"""
function _cos2(U::AbstractMatrix, sigma::AbstractVector, dist2::AbstractVector)
    N, k = size(U)
    @assert length(sigma) == k
    @assert length(dist2) == N
    T = eltype(U)
    out = similar(U)
    @inbounds for j in 1:k, i in 1:N
        d = dist2[i]
        out[i, j] = d > 0 ? T((sigma[j] * U[i, j])^2 / d) : zero(T)
    end
    return out
end

"""
Build the post-SVD output NamedTuple. Centralizing this guarantees identical
field order across `ca` / `sparse_ca` / `bincoo_ca` so the API contract holds.
The first 5 fields preserve the Phase 1/2 positional contract.
"""
function _build_result(U_dim, V_dim, sigma_dim,
                       rowsums, colsums, total,
                       total_inertia, rowdist2, coldist2)

    T = eltype(U_dim)

    inv_sqrt_rowmass = _inv_sqrt_mass(T, rowsums, total)
    inv_sqrt_colmass = _inv_sqrt_mass(T, colsums, total)

    rowmass    = rowsums ./ total                 # Float64 marginal probabilities
    colmass    = colsums ./ total

    rowstd     = _standard_coords(U_dim, inv_sqrt_rowmass)
    colstd     = _standard_coords(V_dim, inv_sqrt_colmass)

    rowcoord   = _principal_coords(rowstd, sigma_dim)
    colcoord   = _principal_coords(colstd, sigma_dim)

    rowcontrib = _contributions(U_dim)
    colcontrib = _contributions(V_dim)

    rowcos2    = _cos2(U_dim, sigma_dim, rowdist2)
    colcos2    = _cos2(V_dim, sigma_dim, coldist2)

    inertia    = sigma_dim .^ 2

    # Phase 4: indices of rows / columns with positive mass. Zero-mass rows
    # and columns flow through the pipeline with all-zero coordinates and
    # zero contributions; users who want to drop them entirely can index
    # into `out.rowcoord[out.valid_rows, :]`.
    valid_rows = findall(>(0), rowsums)
    valid_cols = findall(>(0), colsums)

    # First 5 fields preserve Phase 1/2 positional contract.
    return (rowcoord       = rowcoord,
            colcoord       = colcoord,
            sigma          = sigma_dim,
            inertia        = inertia,
            total_inertia  = total_inertia,
            rowstd         = rowstd,
            colstd         = colstd,
            rowmass        = rowmass,
            colmass        = colmass,
            rowcontrib     = rowcontrib,
            colcontrib     = colcontrib,
            rowcos2        = rowcos2,
            colcos2        = colcos2,
            valid_rows     = valid_rows,
            valid_cols     = valid_cols)
end
