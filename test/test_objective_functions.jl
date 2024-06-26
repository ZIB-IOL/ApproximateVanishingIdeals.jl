using Test
using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals


@testset "Test suite for evaluate_function in L2Loss" begin
  for m in 1:5
    for n in 1:10
      A = rand(m, n)
      b = rand(m)
      x = rand(n)
      lambda = rand() * n
      _, evaluate_function, _ = AVI.L2Loss(A, b, lambda, A' * A, A' * b, b' * b)
            
      @test 1/m * norm(A * x + b, 2)^2 + lambda * norm(x, 2)^2 / 2 ≈ evaluate_function(x)
    end
  end
end;


@testset "Test suite for evaluate_gradient! in L2Loss" begin
  for m in 1:5
    for n in 1:10
      A = rand(m, n)
      b = rand(m)
      x = rand(n)
      lambda = rand() * n
      _, _, evaluate_gradient! = AVI.L2Loss(A, b, lambda, A' * A, A' * b, b' * b)
            
      gradient = 2/m * (A' * A * x + A' * b + m/2 * lambda * x)
      approx_vec = (gradient .≈ evaluate_gradient!(zeros(n), x))
      @test all(approx_vec)
    end
  end
end;
