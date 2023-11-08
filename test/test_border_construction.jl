include("../src/border_construction.jl")
include("../src/auxiliary_functions.jl")

using LinearAlgebra
using Test


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