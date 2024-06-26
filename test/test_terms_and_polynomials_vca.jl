using LinearAlgebra
using ApproximateVanishingIdeals
const AVI = ApproximateVanishingIdeals
using Random
using Test


@testset "Testing suite for apply_V_transformation" begin
    for i in 1:10
        X_train = rand(i, 5)
        X_train_transformed, sets_VCA = AVI.fit_vca(X_train)
        X_transformed, _ = AVI.apply_V_transformation(sets_VCA, X_train)
        X_transformed = abs.(X_transformed)
        @test all(X_train_transformed .- X_transformed .<= 1.0e-10)
    end
end;


@testset "Testing suite for all VCA functions" begin
    X = [   [0.1 0.4]; 
            [0.2 0.5]; 
            [0.3 0.6]   ]
    
    # test construct_SetsVCA
    sets_VCA = AVI.construct_SetsVCA(X)
    @test sets_VCA.X == X
    @test sets_VCA.Cs == []
    @test sets_VCA.Vs == []
    @test sets_VCA.V_coefficient_vectors == []
    @test sets_VCA.Fs[1] == 1/sqrt(size(X, 1)) * ones(Float64, size(X, 1), 1)
    @test sets_VCA.F_coefficient_vectors == [ones(Float64, 1, 1)]
    
    # test construct_border
    border = AVI.construct_border_vca(sets_VCA)
    @test border == X
    
    # test update_F
    F_coeff = [ [1 0.5];
                [0.1 0.2]   ]
    F_eval  = [ [0.5 2];
                [0.3 1];
                [0.4 1]     ]
    AVI.update_F(sets_VCA, F_coeff, F_eval)
    @test sets_VCA.Fs[2] == F_eval
    @test sets_VCA.F_coefficient_vectors[2] == F_coeff
    
    # test update_V
    V_coeff = 0.5 * ones(Float64, 1, 1)
    V_eval = reshape([  [0.1];
                        [0.1];
                        [0.2]   ], 3, 1)
    AVI.update_V(sets_VCA, V_coeff, V_eval)
    @test sets_VCA.Vs[1] == V_eval
    @test sets_VCA.V_coefficient_vectors[1] == V_coeff
    empty_vec = zeros(Float64, 0, 0)
    AVI.update_V(sets_VCA, empty_vec, empty_vec)
    @test sets_VCA.Vs[2] == empty_vec
    @test sets_VCA.V_coefficient_vectors[2] == empty_vec
    V_coeff = reshape([ [0.5];
                        [0.2];
                        [0.0]   ], 3, 1)
    V_eval = reshape([  [-0.3];
                        [0.0];
                        [0.2]   ], 3, 1)
    AVI.update_V(sets_VCA, V_coeff, V_eval)
    @test sets_VCA.Vs[3] == V_eval
    @test sets_VCA.V_coefficient_vectors[3] == V_coeff
    
    # test update_C
    C_evals = reshape([ [0.1];
                        [0.2];
                        [0.3]   ], 3, 1)
    AVI.update_C(sets_VCA, C_evals)
    @test sets_VCA.Cs[1] == C_evals
    
    # test F_to_matrix
    F_coeff = reshape([ [1.0];
                        [-1.0]  ], 2, 1)
    F_eval  = [ [-0.5 0.6];
                [0.3 0.9];
                [0.4 0.8]   ]
    AVI.update_F(sets_VCA, F_coeff, F_eval)
    F_matrix = AVI.F_to_matrix(sets_VCA)
    @test all(F_matrix .- [ [0.57735027 0.5 2.0 -0.5 0.6];
                            [0.57735027 0.3 1.0 0.3 0.9];
                            [0.57735027 0.4 1.0 0.4 0.8]    ] .<= 1.0e-10)
    
    # test V_to_matrix
    V_matrix = AVI.V_to_matrix(sets_VCA)
    @test V_matrix == [ [0.1 -0.3];
                        [0.1 0.0];
                        [0.2 0.2]   ] 
    
    # test evaluate_transformation
    (zeroes, entries, avg_sparsity, number_of_polynomials, number_of_terms, degree) = AVI.evaluate_transformation_vca(sets_VCA)
    @test zeroes == 1
    @test entries == 4
    @test avg_sparsity - 0.25 <= 1.0e-10
    @test number_of_polynomials == 2
    @test number_of_terms == size(AVI.F_to_matrix(sets_VCA), 2)
    @test degree - 2.0 <= 1.0e-10
end;
