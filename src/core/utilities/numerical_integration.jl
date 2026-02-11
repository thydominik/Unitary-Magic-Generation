"""
    numerical_integration

Numerical integration utilities.

Currently this module provides a 2D composite trapezoidal rule integrator for values
sampled on a rectangular grid.

All identifiers and docstrings are ASCII-only.
"""
module numerical_integration

export double_integral_trapz

"""
    double_integral_trapz(f, x, y) -> Float64

Compute a 2D double integral using the composite trapezoidal rule.

The function values are given on a grid defined by vectors `x` and `y`.

Parameters
- `f`: Matrix of size (length(x), length(y)) with f[i,j] = f(x[i], y[j]).
- `x`: Vector of x-axis grid points.
- `y`: Vector of y-axis grid points.

Returns
- Approximate integral value as Float64.

Errors
- Throws DimensionMismatch if the shape of `f` does not match (length(x), length(y)).

Notes
- Grid spacing may be non-uniform.
- The algorithm is O(nx*ny) and allocates O(ny) temporary storage.
"""
function double_integral_trapz(
    f::AbstractMatrix{<:Real},
    x::AbstractVector,
    y::AbstractVector,
)::Float64
    nx, ny = size(f)

    if length(x) != nx || length(y) != ny
        throw(DimensionMismatch(
            "double_integral_trapz: f has size ($nx,$ny) but expected ($(length(x)),$(length(y)))",
        ))
    end

    # Integrate over x for each fixed y.
    integral_x = zeros(Float64, ny)
    @inbounds for j in 1:ny
        acc = 0.0
        for i in 1:(nx - 1)
            dx = Float64(x[i + 1] - x[i])
            acc += 0.5 * dx * (Float64(f[i, j]) + Float64(f[i + 1, j]))
        end
        integral_x[j] = acc
    end

    # Integrate the resulting 1D function over y.
    integral_xy = 0.0
    @inbounds for j in 1:(ny - 1)
        dy = Float64(y[j + 1] - y[j])
        integral_xy += 0.5 * dy * (integral_x[j] + integral_x[j + 1])
    end

    return integral_xy
end

end # module numerical_integration
