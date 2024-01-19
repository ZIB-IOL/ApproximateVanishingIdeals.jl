using Test

include("../src/auxiliary_functions.jl")

matrix = Matrix([[1, 2, 3, 4, 3];; 
                 [0, 1, 2, 3, 2];; 
                 [1, 3, 1, 2, 1];;]')

matrix_sorted = Matrix([[1, 3, 3, 2, 4];; 
                        [0, 2, 2, 1, 3];; 
                        [1, 1, 1, 3, 2];;]')

matrix_unique = Matrix([[1, 3, 2, 4];;
                        [0, 2, 1, 3];;
                        [1, 1, 3, 2];;]')


@testset "Test suite for deg_lex_sort" begin
  matrix_sorted_1, matrix_sorted_2, a = deg_lex_sort(matrix, 1. * matrix)
  @test matrix_sorted_1 == matrix_sorted
  @test matrix_sorted_2 == matrix_sorted
end


@testset "Test suite for get_unique_columns" begin 
  mat1_unique, mat2_unique, _ = get_unique_columns(matrix, matrix)
  @test mat1_unique == matrix_unique
  @test mat2_unique == matrix_unique
end


@testset "Test suite for compute_degree" begin
  @test compute_degree(matrix) == reshape([2, 6, 6, 9, 6], 1, 5)
end
