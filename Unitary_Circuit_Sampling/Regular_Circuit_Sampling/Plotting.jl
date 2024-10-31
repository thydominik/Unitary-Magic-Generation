using Plots
using PlotThemes
# https://github.com/JuliaPlots/PlotThemes.jl

theme(:juno::Symbol;)
f(x) = 10*sin(x)
g(x) = x^2
scatter([f, g], 0, 2*π, title= "F and G", xlabel="x", ylabel="y", linewidth=4)

x = collect(-20:0.1:20);

y1 = x.^2;
y2 = cos.(x);

plot(x, y1, label="x^2", lw=2)
plot(x, y2, label="cos(x)", lw=3)
title!("LinePlot")

#histogram
xlabel!("x")
ylabel!("y")

histogram(y2, alpha = 0.3)

#3D plots