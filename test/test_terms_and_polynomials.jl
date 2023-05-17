using Test

include("../src/terms_and_polynomials.jl")

@testset "Test suite for purge" begin
    matrix_2 = [
        [1, 0, 1],
        [2, 0, 1],
        [3, 1, 1]
    ]
    matrix_1 = [
        [1, 0, 0],
        [1, 1, 1]
    ]
    matrix_2_purged, matrix_2_purged_2, _ = purge(matrix_2, 1.0*matrix_2, matrix_1)
    @test length(matrix_2_purged) == 0
    @test length(matrix_2_purged_2) == 0
    
    matrix_2 = [
        [1, 0, 1],
        [2, 0, 1],
        [3, 1, 1]
    ]
    matrix_1 = [
        [1, 0, 2],
        [1, 2, 1]
    ]
    matrix_2_purged, matrix_2_purged_2, _ = purge(matrix_2, 1.0*matrix_2, matrix_1)
    @test matrix_2_purged == matrix_2_purged_2
    @test matrix_2_purged == [[1, 0, 1], [2, 0, 1], [3, 1, 1]]
end
