begin
    using Plots
    using LinearAlgebra
    using JLD2
    using LaTeXStrings
    using StatsBase
    using MAT
    using ProgressBars
    using StatsBase
    using Random
    using Distributions
end

begin
    function double_integral_trapz(f::Array{<:Real,2}, x::AbstractVector, y::AbstractVector)
        nx, ny = size(f)
        if length(x) != nx || length(y) != ny
            error("Size of f must match lengths of x and y vectors")
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

    function MI(x, y, N_bins)
        m = length(x)

        bins_x      = range(0, 10, N_bins)
        d_bins_x    = bins_x[2] - bins_x[1]

        bins_y      = range(0, 10, N_bins)
        d_bins_y    = bins_y[2] - bins_y[1]

        # Compute joint histogram and probabilities
        joint_hist  = fit(Histogram, (x, y), (bins_x, bins_y))
        pxy = joint_hist.weights ./ (m * d_bins_x * d_bins_y)

        # Marginal probabilities
        px_marginal     = sum(pxy, dims=2)
        h_x             = fit(Histogram, x, bins_x)
        px_individual   = h_x.weights / (m * d_bins_x)

        py_marginal     = sum(pxy, dims=1)
        h_y             = fit(Histogram, y, bins_y)
        py_individual   = h_y.weights / (m * d_bins_y)

        Integrand = zeros(length(px_individual), length(py_individual))
        for i in 1:length(px_individual)
            for j in 1:length(py_individual)
                if pxy[i, j] > 0 
                    Integrand[i, j] = pxy[i, j] * log2(pxy[i, j] / (px_individual[i] * py_individual[j]))
                end
            end
        end
        mi = double_integral_trapz(Integrand, bins_x[1:end-1], bins_y[1:end-1])


        return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, Integrand

    end
end


# STEP 0: Integrate MI numerically and compare that to the analytical solution
NumericalIntegral   = Vector{Float64}()
BinSize             = Vector{Float64}()

rho = 0.25

for m in ProgressBar(1:12)

begin
    mu_x = 0
    mu_y = 0

    sigma_x = 0.5
    sigma_y = 0.1

    rho = rho
    
    p0 = 3
    x = range(-p0 + mu_x, p0 + mu_x, length=2^m)
    y = range(-p0 + mu_y, p0 + mu_y, length=2^m)
end

begin
    Px  = 1/(sqrt(2 * pi) * sigma_x) .* exp.( -0.5 .* ((x .- mu_x) ./ sigma_x).^2 )
    Py  = 1/(sqrt(2 * pi) * sigma_y) .* exp.( -0.5 .* ((y .- mu_y) ./ sigma_y).^2 )

    X = reshape(x, :, 1)  # column vector
    Y = reshape(y, 1, :)  # row vector

    Xc = (X .- mu_x) ./ sigma_x
    Yc = (Y .- mu_y) ./ sigma_y

    Pxy = 1/(2 * pi * sigma_x * sigma_y * sqrt(1 - rho^2)) .* exp.(
        (-1/(2 * (1 - rho^2))) .* (Xc.^2 .+ Yc.^2 .- 2 .* rho .* Xc .* Yc)
    )
end


Integrand = zeros(length(x), length(y))

for i in ProgressBar(1:length(x))
    for j in 1:length(y)
        if Pxy[i, j] > 0
            Integrand[i, j] = Pxy[i, j] * (log2(Pxy[i, j]) -log2((Px[i] * Py[j])))
        end
    end
end

push!(NumericalIntegral, double_integral_trapz(Integrand, x, y))
push!(BinSize, 2*p0/2^m)
end

scatter(BinSize, NumericalIntegral, label = "Numerics")

plot!(xlabel = "Bin size",   xaxis = :log2)
plot!(ylabel = "Analytical MI - Gaussian", yaxis = :lin) 
plot!(box=:on, legend=:topleft, xticks=10)

hline!([-0.5 * log2(1 - rho^2)], label = "Analytics")
vline!([0.01], label="data variance")
savefig("Analytical_MI_Integration_sigma_01_05.png")


# -----------------------------------------------------------------------------


# STEP 1: Gaussian distributions for different ρ - correlation coeff.
#   # How does it converge to the analytics as the sample size M_S increases?
#   # How does it converge to the analytics as M_B or Δ_B?

# No correlation: ρ = 0
begin
    mu_x = 0
    mu_y = 0

    sigma_x = 2
    sigma_y = 2

    rho = .99
    
    p0 = 3
end

CorrelationMatrix   = [sigma_x^2  rho*sigma_x*sigma_y
                        rho*sigma_x*sigma_y  sigma_y^2]
mvnorm              = MvNormal([mu_x, mu_y], CorrelationMatrix)

NumMI = zeros(20, 12)
DataVariance = zeros(20, 2)
Binnings = zeros(20,12, 2)
for m in ProgressBar(1:20)
    M = 2^m
    Data    = rand(mvnorm, M)
    x       = Data[1, :]
    y       = Data[2, :]
    
    DataVariance[m, 1] = var(x)
    DataVariance[m, 2] = var(y)

    for mb in ProgressBar(1:12)

        mesh_x = range(minimum(x), maximum(x), length = 2^mb); dbx = mesh_x[2] - mesh_x[1]
        mesh_y = range(minimum(y), maximum(y), length = 2^mb); dby = mesh_y[2] - mesh_y[1]
        Binnings[m, mb, :] = [dbx, dby]
        hx  = fit(Histogram, x, mesh_x)
        hy  = fit(Histogram, y, mesh_y)
        hxy = fit(Histogram, (x, y), (mesh_x, mesh_y))

        Px  = hx.weights / (M * dbx)
        Py  = hy.weights / (M * dby)
        Pxy = hxy.weights / (M * dbx * dby)

        Integrand = zeros(length(mesh_x) - 1, length(mesh_y) - 1)

        for i in 1:length(mesh_x)-1
            for j in 1:length(mesh_y)-1
                if Pxy[i, j] > 0
                    Integrand[i, j] = Pxy[i, j] * log2(Pxy[i, j] / (Px[i] * Py[j]))
                end
            end
        end

        NumMI[m, mb] = double_integral_trapz(Integrand, mesh_x[1:end-1], mesh_y[1:end-1])

    end
end




NumMI_plot = copy(NumMI)
for i in eachindex(NumMI_plot)
    if isinf(NumMI_plot[i])
        NumMI_plot[i] = NaN
    end
end
heatmap(1:20, 1:12, (NumMI_plot)')
contour!(1:20, 1:12, (NumMI_plot)'; levels=[-0.5 * log2(1 - rho^2)*1.01, -0.5 * log2(1 - rho^2)*0.99], linewidth=2, color=:red)
contour!(1:20, 1:12, (NumMI_plot)'; levels=[-0.5 * log2(1 - rho^2)*1.1 , -0.5 * log2(1 - rho^2)*0.9 ], linewidth=2, color=:blue)
plot!(xlabel = "log2 - Sample size", ylabel = "log2 - Number of Bins", box=:on)
savefig("Gaussian_NumericalMI_BinsVSSamples.png")
#############

plot()
colors = cgrad(:bluesreds, 25, categorical=true)

plot!(yaxis=:log2, ylims=(2^-10, 5))
hline!([-0.5 * log2(1 - rho^2)], color=:black, label="")

plot()
for i in 1:12
    plot!(1:20, NumMI[:, i], label = "log2(M_B) = $(i)")
end
hline!([-0.5 * log2(1 - 0.99^2)], lw = 2, color=:black, label = "analytic solution")
plot!(xlabel = "log2(M_S)", ylabel = "MI", box=:on)
savefig("Gaussian_mutualInformation_test.png")
plot(NumMI)


plot()
for i in 1:20
    plot!(1:12, NumMI[i, :], label = "log2(M_S) = $(i)")
end
hline!([-0.5 * log2(1 - 0.99^2)], lw = 2, color=:black, label = "analytic solution")
plot!(xlabel = "log2(M_B)", ylabel = "MI", box=:on)
savefig("Gaussian_mutualInformation_test_inverted.png")

