# Data for plotting a `LagrangianDescriptorSolution`, shared by the Plots recipe and the
# Makie extension. It returns a `NamedTuple` with the kind of plot and its data:
#
# - `kind = :path`: for scalar initial conditions, the values `z` of the descriptor against
#   the initial conditions `x`;
# - `kind = :heatmap`: for a matrix of initial conditions forming a grid in the components
#   `idxs = (i, j)`, the grid values `x` and `y` of these components, and the values `z` of
#   the descriptor, with `z[k, l]` at `(x[l], y[k])`;
# - `kind = :scatter`: for any other collection of initial conditions, the components `x`
#   and `y` of each initial condition, and the value `z` of the descriptor at each one.
#
# With three components `idxs = (i, j, k)`, for a three-dimensional plot:
#
# - `kind = :surface`: for a matrix of initial conditions, seen as a parametrized surface,
#   the matrices `x`, `y` and `z` of these components, and the matrix `c` of the descriptor;
# - `kind = :scatter3d`: for any other collection of initial conditions, the components `x`,
#   `y` and `z` of each initial condition, and the value `c` of the descriptor at each one.

function _plot_data(lagsol, direction::Symbol, idxs = (1, 2))
    _check_direction(direction, DESCRIPTOR_DIRECTIONS)
    uu0 = lagsol.uu0
    values = lagsol(direction)
    if eltype(uu0) <: Number
        return (kind = :path, x = collect(uu0), y = nothing, z = values)
    end
    first(uu0) isa Union{AbstractVector,Tuple} || throw(
        ArgumentError(
            "Plotting is not implemented for initial conditions of type " *
            "$(typeof(first(uu0))); use numbers, or vectors or tuples of numbers.",
        ),
    )
    n = length(first(uu0))
    if n == 1
        return (kind = :path, x = [only(u0) for u0 in uu0][:], y = nothing, z = values)
    end
    length(idxs) in (2, 3) && all(i -> i isa Integer && 1 ≤ i ≤ n, idxs) || throw(
        ArgumentError(
            "`idxs` must be two or three indices of components of the initial conditions, " *
            "between 1 and $n; got `idxs = $idxs`.",
        ),
    )
    length(idxs) == 3 && return _plot_data_3d(uu0, values, idxs)
    i, j = idxs
    if uu0 isa AbstractMatrix
        Z = reshape(values, size(uu0))
        # `[[x, y] for y in ys, x in xs]`: x varies along the columns and y along the rows
        if _is_grid(uu0, i, j)
            return (kind = :heatmap, x = [uu0[1, c][i] for c in axes(uu0, 2)],
                y = [uu0[r, 1][j] for r in axes(uu0, 1)], z = Z)
        end
        # `[[x, y] for x in xs, y in ys]`: x varies along the rows and y along the columns
        if _is_grid(permutedims(uu0), i, j)
            return (kind = :heatmap, x = [uu0[r, 1][i] for r in axes(uu0, 1)],
                y = [uu0[1, c][j] for c in axes(uu0, 2)], z = permutedims(Z))
        end
    end
    return (kind = :scatter, x = [u0[i] for u0 in uu0][:], y = [u0[j] for u0 in uu0][:], z = values)
end

function _plot_data_3d(uu0, values, (i, j, k))
    if uu0 isa AbstractMatrix && size(uu0, 1) ≥ 2 && size(uu0, 2) ≥ 2
        return (kind = :surface, x = [u0[i] for u0 in uu0], y = [u0[j] for u0 in uu0],
            z = [u0[k] for u0 in uu0], c = reshape(values, size(uu0)))
    end
    return (kind = :scatter3d, x = [u0[i] for u0 in uu0][:], y = [u0[j] for u0 in uu0][:],
        z = [u0[k] for u0 in uu0][:], c = values)
end

# Whether component `i` of the initial conditions only changes along the columns, and
# component `j` only along the rows, with at least two distinct values of each
function _is_grid(uu0::AbstractMatrix, i, j)
    size(uu0, 1) ≥ 2 && size(uu0, 2) ≥ 2 || return false
    for c in axes(uu0, 2), r in axes(uu0, 1)
        uu0[r, c][i] == uu0[1, c][i] || return false
        uu0[r, c][j] == uu0[r, 1][j] || return false
    end
    return uu0[1, 1][i] != uu0[1, 2][i] && uu0[1, 1][j] != uu0[2, 1][j]
end
