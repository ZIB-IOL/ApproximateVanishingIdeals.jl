using Test
using Random
using LinearAlgebra
using FrankWolfe
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals


@testset "Test suite for fit_oavi" begin
  for oracle in ["CG", "Away", "PCG", "Lazy", "BCG", "BPCG"]
    m, n = rand(15:25), rand(4:10)
    X_train = rand(m, n)
    for ihb in ["false", "weak", "full"]
      X_train_transformed, sets = AVI.fit_oavi(X_train; oracle=oracle, inverse_hessian_boost=ihb)
      loss_list = Vector{Float64}([])
      for col in 1:size(sets.G_evaluations, 2)
        cur_col = sets.G_evaluations[:, col]
        loss = 1 / size(cur_col, 1) * (norm(cur_col, 2)^2)
        loss_list = append!(loss_list, loss)
      end
      @test all(loss_list .<= 0.1)
    end
  end
end;  


@testset "Test suite for ABM" begin
  for _ in 1:5
    m, n = rand(15:25), rand(4:10)
    X_train = rand(m, n)
    X_train_transformed, sets = AVI.fit_oavi(X_train; oracle="ABM", psi=0.05)
    loss_list = Vector{Float64}([])
    for col in 1:size(sets.G_evaluations, 2)
      cur_col = sets.G_evaluations[:, col]
      loss = 1 / m * norm(cur_col, 2)^2
      append!(loss_list, loss)
    end
    @test all(loss_list .<= 0.05)
  end
end;


@testset "Test suite for evaluate_oavi (and regularization)" begin
  for oracle in ["CG", "BPCG", "ABM"]
    m, n = rand(15:25), rand(4:10)
    X_tr = rand(m, n)
    X_tr_transformed, sets_tr = AVI.fit_oavi(X_tr; oracle=oracle)
    X_te_transformed, sets_te = AVI.evaluate_oavi(sets_tr, X_tr)

    @test all(X_tr_transformed .- X_te_transformed .<= 1.0e-10)

    if oracle !== "ABM"
        X_train = rand(10, 3)
        X_train_transformed, sets_train = AVI.fit_oavi(X_train; psi=0.01, lambda=0.1, oracle=oracle)
        X_test_transformed, sets_test = AVI.evaluate_oavi(sets_train, X_train)
        @test all(X_train_transformed .- X_test_transformed .<= 1.0e-10)
    end
  end
end;