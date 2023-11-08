using Test
using Random

include("../src/oracle_avi.jl")
include("../src/auxiliary_functions.jl")
include("../src/terms_and_polynomials.jl")
include("../src/objective_functions.jl")
include("../src/auxiliary_functions_avi.jl")


@testset "Test suite for fit" begin
  for oracle in ["CG", "PCG", "BPCG"]
    m, n = rand(15:25), rand(4:10)
    X_train = rand(m, n)
    X_train_transformed, sets = fit(X_train; oracle=oracle)
    loss_list = Vector{Float64}([])
    for col in 1:size(sets.G_evaluations, 2)
      cur_col = sets.G_evaluations[:, col]
      loss = 1 / size(cur_col, 1) * (norm(cur_col, 2)^2)
      loss_list = append!(loss_list, loss)
    end
    @test all(loss_list .<= 0.1)
  end
end;


@testset "Test suite for l1_projection" begin
  vec1 = rand(1:10, 20)
  radius_1 = 2.5
  radius_2 = 3.0
  @test norm(l1_projection(vec1), 1) ≈ 1
  @test norm(l1_projection(vec1; radius=radius_1),1) ≈ radius_1
  @test norm(l1_projection(vec1; radius=radius_2), 1) ≈ radius_2
end
  
