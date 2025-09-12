



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
    using LsqFit
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


begin
    data    = Dict()
    N       = 1:10
    Samples = [1048576 1048576 1048576 1048576 1048576 1048576 442368 1048576 196608 20480]
    for i in 1:10
        data[i] = load("/Volumes/SSD_Szmbthy/PhD_Data/Unitary_Circuits_Magic_and_Entanglement/Regular/RegularUnitaryCircuitMagicSampled_N_$(N[i])_Samples_$(Samples[i]).jld2")
    end
end

# N = 2 case

MutualInformation = zeros(21, 13)
DataVarianceM    = Vector{Float64}()
DataVarianceS   = Vector{Float64}()

for ms in ProgressBar(21)
    M = data[2]["Magic"][1:2^ms] .* log(3/2) ./ log2(3/2);   SampleSizeM = length(M)
    push!(DataVarianceM, var(M))
    S = data[2]["Svn"][1:2^ms];      SampleSizeS = length(S)
    push!(DataVarianceS, var(S))
    for m in 8

        binsM = range(0, log2(3/2), length=2^m);    dbM = binsM[2] - binsM[1]
        binsS = range(0, 1, length=2^m);            dbS = binsS[2] - binsS[1]

        hx  = fit(Histogram, M, binsM)
        hy  = fit(Histogram, S, binsS)
        hxy = fit(Histogram, (M, S), (binsM, binsS))

        Px = hx.weights / (SampleSizeM * dbM)
        Py = hy.weights / (SampleSizeS * dbS)
        Pxy = hxy.weights / (SampleSizeM * dbM * dbS)

        global Integrand = zeros(length(binsM) - 1, length(binsS) - 1)
        for i in 1:2^m - 1
            for j in 1:2^m - 1
                if Pxy[i, j] > 0
                    Integrand[i, j] = Pxy[i, j] * log2(Pxy[i, j] / (Px[i] * Py[j]))
                end
            end
        end

        MutualInformation[ms, m] = double_integral_trapz(Integrand, binsM[1:end-1], binsS[1:end-1])

    end
end

heatmap(1:21, 1:13, MutualInformation')
contour!(1:21, 1:13, (MutualInformation)';levels = [2^-5, 2^-4, 2^-3, 2^-2], linewidth=2, color=:red)

binNumbers = 2 .^ (1:13)
binSizeM = log2(3/2) ./ binNumbers
binSizeS = 1 ./ binNumbers
scatter(binSizeM)
hline!([DataVarianceM[21]], yaxis=:log2)



N   = 5
NMS = 15
NMB = 10

Magic = data[N]["Magic"];   SampleSizeM = length(Magic)
Svn = data[N]["Svn"];     SampleSizeS = length(Svn)

Minimum_M = minimum(Magic)
Maximum_M = maximum(Magic)

Minimum_S = minimum(Svn)
Maximum_S = maximum(Svn)

MutualInformation   = zeros(NMS, NMB)
DataVarianceM       = Vector{Float64}()
DataVarianceS       = Vector{Float64}()

for ms in ProgressBar(5:NMS)  
    if ms == NMS
        M = Magic[1:end]; SampleSizeM = length(M)
        S = Svn[1:end]; SampleSizeS = length(S)
    else
        M = Magic[1:2^ms]; SampleSizeM = length(M)
        S = Svn[1:2^ms]; SampleSizeS = length(S)
    end

    

    push!(DataVarianceM, var(M))
    push!(DataVarianceS, var(S))

    for m in 5:NMB

        binsM = range(Minimum_M, Maximum_M, length=2^m);    dbM = binsM[2] - binsM[1]
        binsS = range(Minimum_S, Maximum_S, length=2^m);    dbS = binsS[2] - binsS[1]

        hx  = fit(Histogram, M, binsM)
        hy  = fit(Histogram, S, binsS)
        hxy = fit(Histogram, (M, S), (binsM, binsS))

        Px = hx.weights / (sum(hx.weights) * dbM)
        Py = hy.weights / (sum(hy.weights) * dbS)
        Pxy = hxy.weights / (sum(hxy.weights) * dbM * dbS)

        global Integrand = zeros(length(binsM) - 1, length(binsS) - 1)
        for i in 1:(2^m - 1)
            for j in 1:(2^m - 1)
                if Pxy[i, j] > 0
                    Integrand[i, j] = Pxy[i, j] * log2(Pxy[i, j] / (Px[i] * Py[j]))
                end
            end
        end

        MutualInformation[ms, m] = double_integral_trapz(Integrand, binsM[1:end-1], binsS[1:end-1])

    end
end

heatmap(MutualInformation')
contour!((MutualInformation)';levels = [2^-8, 2^-4, 2^-3, 2^-2], linewidth=2, color=:red)


plot(MutualInformation[:, 8])
plot!(MutualInformation[:, 7])
plot!(MutualInformation[:, 5])





# Select the last few points for fitting (e.g., last 6 points)
xdata = Float64.(collect(size(MutualInformation, 1)-5:size(MutualInformation, 1))); xdata[end] = log2(20480) 
ydata = MutualInformation[end-5:end, 6]

# Exponential model: y = a * exp(b * x) + c
exp_model(x, p) = p[1] .* exp.(x .* (p[2])) .+ p[3]  # Ensure offset is positive

# Initial guess for parameters [a, b, c]
p0 = [ydata[end], -0.1, 0.00]  # Set initial offset positive

fitt = LsqFit.curve_fit(exp_model, xdata, ydata, p0)
fitted_params = fitt.param

#println("Fitted parameters: a=$(fitted_params[1]), b=$(fitted_params[2]), c=$(fitted_params[3]), d=$(fitted_params[4])")
println("Fitted parameters: a=$(fitted_params[1]), b=$(fitted_params[2]), c=$(fitted_params[3])")
stderror(fitt)
# Plot original data and fit
plot(xdata, ydata, label="Data", marker=:o)
plot!(xdata, exp_model(xdata, fitted_params), label="Exponential Fit", lw=2)

Fitting_And_Error = zeros(4, 2)
for p in 5:8
    ydata = MutualInformation[end-5:end, p]

    p0 = [ydata[end], -0.1, 0.00]

    fitt = LsqFit.curve_fit(exp_model, xdata, ydata, p0)
    fitted_params = fitt.param

    Fitting_And_Error[p - 4, 1] = fitted_params[3]
    Fitting_And_Error[p - 4, 2] = stderror(fitt)[3]
end

@save "MI_fitting_exponential_N$(N)_BinningSizeDependence.jld2" Fitting_And_Error



Params = zeros(10, 4, 2)
for n in 2:10
    global data = load("MI_fitting_exponential_N$(n)_BinningSizeDependence.jld2")["Fitting_And_Error"]
    Params[n, :, :] = data
end 

plot(2:10, abs.(Params[2:end, 1, 1]), yerror=Params[2:end, 1, 2], xlims=(1, 11), yaxis=:log2)
plot!(2:10, abs.(Params[2:end, 2, 1]), xlims=(2, 10), yaxis=:log2)
plot!(2:10, abs.(Params[2:end, 3, 1]), yerror=Params[2:end, 3, 2],xlims=(2, 10), yaxis=:log2)
plot!(2:10, abs.(Params[2:end, 4, 1]), yerror=Params[2:end, 4, 2],xlims=(2, 10), yaxis=:log2)


Data = Params[:, 3, :]
@save "MB_8_Fitting_P3_with error.jld2" Data

A = load("MB_7_Fitting_P3_with error.jld2")