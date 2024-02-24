using Test
using Random
using LinearAlgebra
using FrankWolfe

include("../src/oracle_avi.jl")
include("../src/auxiliary_functions.jl")
include("../src/terms_and_polynomials.jl")
include("../src/objective_functions.jl")
include("../src/auxiliary_functions_avi.jl")
include("../src/oracle_constructors.jl")
include("../src/border_construction.jl")

@testset "Test suite for fit_oavi" begin
  for oracle in ["CG", "BCG", "BPCG"]
    m, n = rand(15:25), rand(4:10)
    X_train = rand(m, n)
    for ihb in ["false", "weak", "full"]
      X_train_transformed, sets = fit_oavi(X_train; oracle=oracle, inverse_hessian_boost=ihb)
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
    X_train_transformed, sets = fit_oavi(X_train; oracle="ABM", psi=0.05)
    loss_list = Vector{Float64}([])
    for col in 1:size(sets.G_evaluations, 2)
      cur_col = sets.G_evaluations[:, col]
      loss = 1 / m * norm(cur_col, 2)^2
      append!(loss_list, loss)
    end
    @test all(loss_list .<= 0.05)
  end
end;
