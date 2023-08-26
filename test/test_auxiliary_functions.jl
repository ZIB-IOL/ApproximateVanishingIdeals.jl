using Test

include("../src/auxiliary_functions.jl")


@testset "Test suite for vecvec_to_mat" begin
  vecvec1 = [[1, 2, 3], [4, 5, 6]]
  vecvec2 = [[1, 2], [3, 4], [5, 6]]

  row_1 = reshape([1, 2, 3], 1, 3)
  row_2 = reshape([4, 5, 6], 1, 3)
  
  mat1 = vcat(row_1, row_2)
  mat2 = reshape([1, 2, 3, 4, 5, 6], 2, 3)

  @test vecvec_to_mat(vecvec1) == mat1
  @test vecvec_to_mat(vecvec2, arr_is_col=1) == mat2
end


matrix = vecvec_to_mat([[1, 2, 3, 4, 3], [0, 1, 2, 3, 2], [1, 3, 1, 2, 1]])
matrix_sorted = vecvec_to_mat([[1, 3, 3, 2, 4], [0, 2, 2, 1, 3], [1, 1, 1, 3, 2]])
matrix_unique = vecvec_to_mat([[1, 3, 2, 4], [0, 2, 1, 3], [1, 1, 3, 2]])


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
