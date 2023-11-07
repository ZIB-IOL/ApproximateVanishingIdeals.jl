using Test

include("../src/auxiliary_functions.jl")
include("../src/terms_and_polynomials.jl")


@testset "Test suite for purge" begin
  matrix_2 = vecvec_to_mat([[1, 2, 3], 
                            [0, 0, 1],
                            [1, 1, 1]])
    
  matrix_1 = vecvec_to_mat([[1, 1],
                            [0, 1],
                            [0, 1]])
    
  matrix_2_purged, matrix_2_purged_2, _ = purge(matrix_2, 1. * matrix_2, matrix_1)
    
  @test size(matrix_2_purged, 2) == 0
  @test size(matrix_2_purged_2, 2) == 0
    
  matrix_2 = vecvec_to_mat([[1, 2, 3],
                            [0, 0, 1],
                            [1, 1, 1]])
    
  matrix_1 = vecvec_to_mat([[1, 1],
                            [0, 2],
                            [2, 1]])
    
  matrix_2_purged, matrix_2_purged_2, _ = purge(matrix_2, 1. * matrix_2, matrix_1)
    
  @test matrix_2_purged == matrix_2_purged_2
  @test matrix_2_purged == vecvec_to_mat([[1, 2, 3],
                                          [0, 0, 1],
                                          [1, 1, 1]])    
end


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
    
  @test G_coefficient_vectors == vecvec_to_mat([[1., 1., 1., 1.],
                                                [2., 0., 0., 0.],
                                                [0., 2., 0., 0.],
                                                [0., 0., 2., 0.],
                                                [0., 0., 0., 2.],
                                                [0., 0., 0., 3.]])
end


@testset "Test suite for apply_G_transformation" begin
    for i in 1:5
        X_train = rand(2*i, 3)
        X_train_transformed, sets_train = fit(X_train)
        
        if X_train_transformed != nothing
            X_test_transformed, sets_test = apply_G_transformation(sets_train, X_train)
            @test X_test_transformed != nothing
            
            @test all(abs.(X_test_transformed) .- X_train_transformed .<= 1.0e-10)
        end
    end
end;
