using Test
using FrankWolfe
using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals


@testset "Test suite for apply_G_transformation" begin
    for i in 1:5
        X_train = rand(2*i, 3)
        X_train_transformed, sets_train = AVI.fit_oavi(X_train; psi=0.01)
        
        if X_train_transformed !== nothing
            X_test_transformed, sets_test = AVI.apply_G_transformation(sets_train, X_train)
            @test X_test_transformed !== nothing            
            @test all(abs.(X_test_transformed) .- X_train_transformed .<= 1.0e-10)
        end
    end
end;


@testset "Test suite for evaluate_transformation_oavi" begin
    X_train = [ [1 2];
                [3 4];
                [5 6]   ]
    
    sets_avi = AVI.construct_SetsOandG(X_train)

    sets_avi.G_coefficient_vectors = [nothing, nothing, [   [1 0];
                                                            [2 1];
                                                            [3 0];
                                                            [0 3]   ], nothing, [   [1 1 1 1];
                                                                                    [2 0 0 0];
                                                                                    [0 1 0 0];
                                                                                    [0 0 2 0];
                                                                                    [0 0 1 0];
                                                                                    [0 0 0 1]   ]]

    (total_number_of_zeros, total_number_of_entries, avg_sparsity, number_of_polynomials, number_of_terms, degree
        ) = AVI.evaluate_transformation_oavi(sets_avi)

    @test total_number_of_zeros == 2
    @test total_number_of_entries == 16
    @test abs(avg_sparsity - 0.16666666666666) < 1.0e-8
    @test number_of_terms == size(sets_avi.O_evaluations, 2)
    @test number_of_polynomials == 6
    @test abs(degree - 3.333333333333) < 1.0e-8         
end
