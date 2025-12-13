"""
    NumericalIntegration

Module for numerical integration utilities, primarily for 2D trapezoidal integration.
"""
module NumericalIntegration

export double_integral_trapz

"""
    double_integral_trapz(f::Array{<:Real,2}, x::AbstractVector, y::AbstractVector)::Float64

Computes a 2D double integral using the trapezoidal rule.

# Arguments
- `f::Array{<:Real,2}`: 2D array of function values
- `x::AbstractVector`: x-axis grid points
- `y::AbstractVector`: y-axis grid points

# Returns
- `Float64`: The computed integral value

# Throws
- `error` if dimensions don't match

# Example
```julia
f = randn(100, 100)
x = range(0, 1, 100)
y = range(0, 1, 100)
result = double_integral_trapz(f, x, y)
```
"""
function double_integral_trapz(f::Array{<:Real,2}, x::AbstractVector, y::AbstractVector)::Float64
    nx, ny = size(f)
    if length(x) != nx || length(y) != ny
        error("Size of f ($(size(f))) must match lengths of x ($(length(x))) and y ($(length(y))) vectors")
    end

    # Trapezoidal integral over x for each fixed y
    integral_x = zeros(ny)
    for j in 1:ny
        for i in 1:(nx-1)
            dx = x[i+1] - x[i]
            integral_x[j] += 0.5 * dx * (f[i,j] + f[i+1,j])
        end
    end

    # Trapezoidal integral over y of the results
    integral_xy = 0.0
    for j in 1:(ny-1)
        dy = y[j+1] - y[j]
        integral_xy += 0.5 * dy * (integral_x[j] + integral_x[j+1])
    end

    return integral_xy
end

end  # module NumericalIntegration
