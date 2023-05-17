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


@testset "Test suite for monomial_evaluation_set" begin
    m1 = [4, 1, 2]
    X1 = [[1.0, 2.0, 3.0], [1.0, 0.5, 2.0]]
    coeff = 4.5
    @test monomial_evaluation_set(m1, X1) == [18.0, 2.0]
    @test monomial_evaluation_set(m1, X1, coeff) == [81.0, 9.0]
end


@testset "Test suite for monomial_set_evaluation_set" begin
    M1 = [[4, 1, 2], [3, 0, 2], [0, 2, 1]]
    X1 = [[1.0, 2.0, 3.0], [1.0, 0.5, 2.0]]
    coeffs = [2.0, 1.0, 5.0]
    @test monomial_set_evaluation_set(M1, X1) == [[18.0, 2.0], [9.0, 4.0], [12.0, 0.5]]
    @test monomial_set_evaluation_set(M1, X1, coeffs) == [[36.0, 4.0], [9.0, 4.0], [60.0, 2.5]]
end
