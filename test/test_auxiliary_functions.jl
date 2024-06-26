using Test
using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals


matrix = Matrix([[1 2 3 4 3]; 
                 [0 1 2 3 2]; 
                 [1 3 1 2 1];])

matrix_sorted = Matrix([[1 3 3 2 4]; 
                        [0 2 2 1 3]; 
                        [1 1 1 3 2];])

matrix_unique = Matrix([[1 3 2 4];
                        [0 2 1 3];
                        [1 1 3 2];])


@testset "Test suite for deg_lex_sort" begin
  matrix_sorted_1, matrix_sorted_2, _ = AVI.deg_lex_sort(matrix, 1. * matrix)
  @test matrix_sorted_1 == matrix_sorted
  @test matrix_sorted_2 == matrix_sorted
end


@testset "Test suite for get_unique_columns" begin 
  mat1_unique, mat2_unique, unique_inds = AVI.get_unique_columns(matrix, matrix)
  @test mat1_unique == matrix_unique
  @test mat2_unique == matrix_unique
  @test unique_inds == [1, 2, 4, 5]
  
  mat1_unique, mat2_unique, _ = AVI.get_unique_columns(matrix)
  @test mat1_unique == matrix_unique
  @test mat2_unique == zeros(Float64, 0, 0)
end


@testset "Test suite for compute_degree" begin
  @test AVI.compute_degree(matrix) == [2 6 6 9 6]
end


matrix_non_zero = Matrix([[1 0 3 0 0 0]; 
                          [0 1 2 0 0 0]; 
                          [1 3 1 2 0 0];
                          [1 0 0 0 3 0]])

@testset "Test suite for finding non-zero entries" begin
  first_ids = AVI.find_first_non_zero_entries(matrix_non_zero)
  @test first_ids == [1, 2, 1, 3, 4, 1]

  last_ids = AVI.find_last_non_zero_entries(matrix_non_zero)
  @test last_ids == [4, 3, 3, 3, 4, 4]
end