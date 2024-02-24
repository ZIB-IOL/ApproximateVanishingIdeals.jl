using Test
using FrankWolfe

include("../src/auxiliary_functions.jl")
include("../src/terms_and_polynomials.jl")
include("../src/oracle_constructors.jl") 
include("../src/border_construction.jl")
include("../src/objective_functions.jl")
include("../src/auxiliary_functions_avi.jl")
include("../src/oracle_avi.jl")

@testset "Test suite for update_coefficient_vectors" begin
  G_coefficient_vectors = reshape([1, 2, 0], 3, 1)
    
  vec1 = reshape([1, 2], 2, 1)
    
  G_coefficient_vectors = update_coefficient_vectors(G_coefficient_vectors, vec1)
  G_coefficient_vectors = vcat(G_coefficient_vectors, zeros(1, size(G_coefficient_vectors, 2)))
  G_coefficient_vectors = update_coefficient_vectors(G_coefficient_vectors, vec1)
  G_coefficient_vectors = vcat(G_coefficient_vectors, zeros(1, size(G_coefficient_vectors, 2)))
  G_coefficient_vectors = vcat(G_coefficient_vectors, zeros(1, size(G_coefficient_vectors, 2)))
    
  vec1 = reshape([1, 2, 3], 3, 1)
    
  G_coefficient_vectors = update_coefficient_vectors(G_coefficient_vectors, vec1)
    
  @test G_coefficient_vectors == Matrix([[1. 1. 1. 1.];
                                         [2. 0. 0. 0.];
                                         [0. 2. 0. 0.];
                                         [0. 0. 2. 0.];
                                         [0. 0. 0. 2.];
                                         [0. 0. 0. 3.];])
end


@testset "Test suite for l1_projection" begin
    vec1 = rand(1:10, 20)
    radius_1 = 2.5
    radius_2 = 3.0
    @test norm(l1_projection(vec1), 1) ≈ 1
    @test norm(l1_projection(vec1; radius=radius_1), 1) ≈ radius_1
    @test norm(l1_projection(vec1; radius=radius_2), 1) ≈ radius_2

    vec2 = l1_projection(vec1)
    @test vec2 ≈ l1_projection(vec2)
    @test norm(l1_projection(vec2), 1) ≈ 1
    @test norm(vec2) ≈ norm(l1_projection(vec2))

  end