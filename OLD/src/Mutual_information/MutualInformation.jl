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


function MIn(x, y, N_bins)
    m = length(x)

    bins_x      = range(minimum(x), maximum(x), N_bins)
    d_bins_x    = bins_x[2] - bins_x[1]

    bins_y      = range(minimum(y), maximum(y), N_bins)
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
                Integrand[i, j] = pxy[i, j] / m * log2(pxy[i, j] * m / (px_individual[i] * py_individual[j]))
            end
        end
    end
    mi = double_integral_trapz(Integrand, bins_x[1:end-1], bins_y[1:end-1])
    

    return mi, px_individual, px_marginal, py_individual, py_marginal, pxy, Integrand

end


end

data    = Dict()
N       = 1:10
Samples = [1048576 1048576 1048576 1048576 1048576 1048576 442368 1048576 196608 20480]
for i in 1:10
    data[i] = load("/Volumes/SSD_Szmbthy/PhD_Data/Unitary_Circuits_Magic_and_Entanglement/Regular/RegularUnitaryCircuitMagicSampled_N_$(N[i])_Samples_$(Samples[i]).jld2")
end




for N in ProgressBar(2:10)
    x = data[N]["Magic"]
    y = data[N]["Svn"]

    MIM     = Vector{Float64}()
    MS      = Vector{Float64}()

    for m in ProgressBar(0:14)
        selected_indices = randperm(length(x))[1:2^m]
        MutualInformation, ~ = MI(x[selected_indices], y[selected_indices], 2^12)
        push!(MIM, MutualInformation)
        push!(MS, 2^m)
    end

    result = Dict("MIM" => MIM, "MS" => MS)
    @save "MIn_N$(N)_sampleSizeDependence.jld2" result
end

plot()
scatter!(log2.(MS), MIM, label = "N = $N")
plot!(xlabel = "Sample size - log2(M)")
plot!(ylabel = "Mutual Information")
plot!(box = :on)
savefig("MI_N$(N)_sampleSizeDependence.png")



all_MIM = []
all_MS = []
labels = []

for n in 2:10
    fname = "MI_N$(n)_sampleSizeDependence.jld2"
    if isfile(fname)
        result = load(fname)
        push!(all_MIM,  result["result"]["MIM"])
        push!(all_MS,   result["result"]["MS"])
        push!(labels, "N = $n")
    end
end

plot()

colors = [RGB(0, 0, 1) * (1 - (i-1)/(length(all_MIM)-1)) + RGB(1, 0, 0) * ((i-1)/(length(all_MIM)-1)) for i in 1:length(all_MIM)]

for i in 1:length(all_MIM)
    ms_vals = all_MS[i]
    mim_vals = all_MIM[i]
    unique_ms = unique(ms_vals)
    global avg_mim = [mean(mim_vals[ms_vals .== ms]) for ms in unique_ms]
    scatter!((unique_ms), avg_mim, label = labels[i], color = colors[i])
    println(avg_mim)

end

plot!(xlabel = "Sample size - log2(M)")
plot!(ylabel = "Mutual Information")
plot!(box = :on)
lambda = .2;
x = range(0, 22, 1000)
func = 2 .* x .* exp.(-lambda .* x)
plot!(xaxis=:log2)
savefig("MI_AllN_sampleSizeDependence.png")
    
plot()
mim_at_2_14 = [mean(all_MIM[i][all_MS[i] .== 2^20]) for i in 1:length(all_MIM)]
ns = 2:10
scatter!(ns, mim_at_2_14, label = "MIM at MS = 2^20", color = :blue)
plot!(xlabel = "N")
plot!(ylabel = "Mutual Information")
plot!(box = :on, yaxis=:lin)
savefig("MI_AllN_MS_2^20.png")



N = 5
x = data[N]["Magic"]
y = data[N]["Svn"]

MIM = Vector{Float64}()
MIMn = Vector{Float64}()
MBS = Vector{Float64}()
for i in ProgressBar(1:14)
    
    selected_indices = randperm(length(x))
    MInfn, ~ = MIn(x[selected_indices], y[selected_indices], 2^i)
    MInf, ~ = MI(x[selected_indices], y[selected_indices], 2^i)
    push!(MIMn, MInfn)
    push!(MIM, MInf)
end


result = Dict("MIM" => MIM, "MBS" => MBS)
@save "MI_N$(N)_BinningSizeDependence.jld2" result

scatter(MIM, label = "N = 5 with M = 2^10")
scatter(MIMn)
plot!(xlabel = "Number of bins - log2(M_B)", ylabel = "MI")
vline!([log2(10000)], label = "used number of bins in the calculations: 10000", color = :red, linestyle = :dash)
savefig("MI_N$(N)_BinNumberDependence_MS_2^10.png")




## Quick Gaussian test

μ = 5           # mean
σ = 0.5           # standard deviation
n = 10^4          # number of samples

Gauss_x = randn(n) .* σ .+ μ
histogram(Gauss_x, bins=range(0, 10, length=100))

Gauss_y_dep = Gauss_x
Gauss_y_indep = randn(n) .* σ .+ μ
histogram!(Gauss_y_indep, bins=range(0, 10, length=100))

mi, px_individual, px_marginal, py_individual, py_marginal, pxy = MI(Gauss_x, Gauss_y_indep, 10000)
mi, px_individual, px_marginal, py_individual, py_marginal, pxy = MI(Gauss_x, Gauss_y_dep, 100)






## Histogram normalization test
Ms = 10^4
x = randn(Ms)
histogram(x)
bins = range(-5, 5, length = 1000)

# probability of x: P(x)
hx = fit(Histogram, x, bins)
P_x = hx.weights ./ Ms
sum(P_x)

# Probability density function
dx = bins[2] - bins[1]
Rho_x = hx.weights ./ ((Ms) * dx)
sum(Rho_x)
Integral = 0.0
for i in 1:length(Rho_x)-1
    Integral += dx * (Rho_x[i] + Rho_x[i + 1])/2
end 
Integral


# Probability density function of P(x, y)
Ms = 10^5
x = randn(Ms) .+ 4
y = randn(Ms) .+ 5
bins = range(0, 10, length=1000)
db = bins[2] - bins[1]

hxy = fit(Histogram, (x, y), (bins, bins))

Pxy = hxy.weights ./ Ms
sum(Pxy)

Rxy = hxy.weights ./ (Ms * db^2)


double_integral_trapz(Rxy, bins[1:end-1], bins[1:end-1])

mi, ~ = MI(x, x, 10000)
0.5 * log(2*pi*exp(1))


