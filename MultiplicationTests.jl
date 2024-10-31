using LinearAlgebra, BenchmarkTools

U = rand(Float64, 1024, 1024) + rand(Float64, 1024, 1024);
psi = ones(2^10);

NPsi = U * psi;