using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals
using Test


@testset "Test suite for purge" begin
  matrix_2 = Matrix([[1 2 3]; 
                     [0 0 1];
                     [1 1 1];])
    
  matrix_1 = Matrix([[1 1];
                     [0 1];
                     [0 1];])
    
  matrix_2_purged, matrix_2_purged_2, _ = AVI.purge(matrix_2, 1. * matrix_2, matrix_1)
    
  @test size(matrix_2_purged, 2) == 0
  @test size(matrix_2_purged_2, 2) == 0
    
  matrix_2 = Matrix([[1 2 3];
                     [0 0 1];
                     [1 1 1];])
    
  matrix_1 = Matrix([[1 1];
                     [0 2];
                     [2 1];])
    
  matrix_2_purged, matrix_2_purged_2, _ = AVI.purge(matrix_2, 1. * matrix_2, matrix_1)
    
  @test matrix_2_purged == matrix_2_purged_2
  @test matrix_2_purged == Matrix([[1 2 3];
                                   [0 0 1];
                                   [1 1 1];])    
end


@testset "Test suite for construct_border" begin
  degree_1_terms = 1 * Matrix(I, 3, 3)

  terms = hcat(degree_1_terms, [1, 0, 1], [0, 1, 1])

  # only purging term is [0, 1, 1] = x2x3
  purging_terms = hcat(zeros(Int64, 3, 0), terms[:, 5])

  # how the raw border looks
  raw_border = [  [2  1  1  0  1  1  0  0  0  2  1  1  0  1  0];
                  [0  1  1  2  0  0  1  1  0  0  1  1  2  0  1];
                  [0  0  0  0  1  1  1  1  2  1  1  1  1  2  2]   ]

  # duplicate indices: 3, 6, 8, 12; purged indices: 7, 8, 11, 12, 13, 15
  unique_non_purging_indices = [1, 2, 4, 5, 9, 10, 14]

  terms_raw, _, non_purging_indices, _ = AVI.construct_border(terms, 1. * terms, zeros(Float64, 0, 0), degree_1_terms, 1. * degree_1_terms, purging_terms)

  @test terms_raw[:, non_purging_indices] == raw_border[:, unique_non_purging_indices]
end