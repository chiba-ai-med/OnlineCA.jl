"""
Supplementary point projection into a trained CA space.

For an active row `i` with row profile `p_i = P_i / r_i`, the principal row
coordinate is

    F_ik = Σ_j (p_ij - c_j) * V_jk / sqrt(c_j) = (p_i - c)' * colstd[:, k]

(The derivation factors S V = U Σ through the D_r^{-1/2} normalization;
see Greenacre 2017 §12.) The same formula applied to a supplementary
(passive) row `x*` projects it into the existing CA space without
re-fitting. By symmetry, supplementary columns project through `rowstd`.

These functions are deliberately format-agnostic: they accept the result
NamedTuple and a small in-memory matrix of new points. An OOC variant
that streams `X_new` from disk is a straightforward future extension.
"""

"""
    project_rows(result, X_new)

Project supplementary rows into the CA space spanned by `result`. Each row
of `X_new` is one supplementary point with the same column count as the
training data. Returns an `n_sup × dim` matrix of principal row
coordinates.
"""
function project_rows(result, X_new::AbstractMatrix)
    M = length(result.colmass)
    size(X_new, 2) == M ||
        throw(DimensionMismatch("X_new has $(size(X_new, 2)) columns; expected $M"))
    rowtotals = sum(X_new, dims=2)
    any(iszero, rowtotals) &&
        throw(ArgumentError("supplementary rows must have a positive total"))
    Rprof = X_new ./ rowtotals                          # n_sup × M
    return (Rprof .- result.colmass') * result.colstd  # n_sup × dim
end

"""
    project_columns(result, Y_new)

Project supplementary columns. Each column of `Y_new` is one supplementary
point with the same row count as the training data. Returns an
`n_sup × dim` matrix of principal column coordinates.
"""
function project_columns(result, Y_new::AbstractMatrix)
    N = length(result.rowmass)
    size(Y_new, 1) == N ||
        throw(DimensionMismatch("Y_new has $(size(Y_new, 1)) rows; expected $N"))
    coltotals = sum(Y_new, dims=1)
    any(iszero, coltotals) &&
        throw(ArgumentError("supplementary columns must have a positive total"))
    Cprof = Y_new ./ coltotals                          # N × n_sup
    return (Cprof .- result.rowmass)' * result.rowstd  # n_sup × dim
end
