include("../src/terms_and_polynomials_vca.jl");
include("../src/vca.jl");
include("../src/auxiliary_functions.jl");

using LinearAlgebra
using Random
using Test


@testset "Testing suite for VCA fit" begin
    psi = 0.1
    max_degree = 10
    for i in 1:5
        X_train = rand(i, 3)
        X_train_transformed, sets_VCA = fit_vca(X_train; psi=psi, max_degree=max_degree)
        for j in 1:size(X_train_transformed, 2)
            @test 1/i * norm(X_train_transformed[:, j], 2)^2 <= psi
        end
    end
end;


@testset "Testing suite for find_range_null_vca" begin
    for run in 1:10
        F = rand(5, run)
        C = rand(5, 2*run)
        psi = 0.1
        
        V_coeffs, V_evals, F_coeffs, F_evals = find_range_null_vca(F, C, psi)
        
        for i in 1:size(F_coeffs, 2)
            val = abs.(hcat(F, C) * F_coeffs[:, i] - F_evals[:, i])
            @test all(val .<= 10e-10)
            @test any(abs.(F_evals[:, i]) .> psi) 
        end
        
        for i in 1:size(V_coeffs, 2)
            val = abs.(hcat(F, C) * V_coeffs[:, i] - V_evals[:, i])
            @test all(val .<= 10e-10)
            @test all(1/5 .* norm(V_evals[:, i])^2 .<= psi)
        end
        
    end
end;


@testset "Testing suite for evaluate_vca" begin 
    m, n = rand(1:50), rand(1:50)
    X_train = rand(m, n)
    X_train_transformed, sets_train = fit_vca(X_train)
    
    if X_train_transformed != nothing
        X_test_transformed, sets_test = evaluate_vca(sets_train, X_train)
        @test X_test_transformed != nothing
        @test all(X_train_transformed .- X_test_transformed .<= 1.0e-10)
        
        F_train_coeffs = sets_train.F_coefficient_vectors
        F_test_coeffs = sets_test.F_coefficient_vectors
        @test length(F_train_coeffs) == length(F_test_coeffs)
        
        for idx in 1:length(F_train_coeffs)
            F_train_coeffs_cp = F_train_coeffs[idx]
            F_test_coeffs_cp = F_test_coeffs[idx]
            @test all(abs.(F_train_coeffs_cp .- F_test_coeffs_cp) .<= 1.0e-10)
        end
        
        V_train_coeffs = sets_train.V_coefficient_vectors
        V_test_coeffs = sets_test.V_coefficient_vectors
        @test length(V_train_coeffs) == length(V_test_coeffs)
        
        for idx in 1:length(V_train_coeffs)
            if size(V_train_coeffs[idx], 2) == 0
                @test size(V_test_coeffs[idx], 2) == 0
            else
                V_train_coeffs_cp = V_train_coeffs[idx]
                V_test_coeffs_cp = V_test_coeffs[idx]
                @test all(abs.(V_train_coeffs_cp .- V_test_coeffs_cp) .<= 1.0e-10)
            end
        end
    end
end;
