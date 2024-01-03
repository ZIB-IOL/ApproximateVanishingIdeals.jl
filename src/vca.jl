"""
This function creates and applies a VCA transformation fitted to X.

# Arguments
- 'X::Matrix{Float64}': data, stored row-wise
- 'psi::Float64': vanishing parameter
- 'max_degree::Int64': maximum degree to consider

# Returns
- 'X_train_transformed::Matrix{Float64}': X transformed according to transformation found by VCA
- 'sets_vca::SetsVCA': instance of SetsVCA containing relevant sets for VCA
"""
function fit_vca(X::Matrix{Float64};
    psi::Float64=0.1,
    max_degree::Int64=10)
        
    sets_vca = construct_SetsVCA(X)
    degree = 1
    while degree < max_degree
        border = construct_border(sets_vca; degree=degree)
        
        update_C(sets_vca, border)
        
        V_coefficients, V_evaluations, F_coefficients, F_evaluations = find_range_null_vca(
            F_to_matrix(sets_vca), sets_vca.Cs[end], psi)
        
        update_V(sets_vca, V_coefficients, V_evaluations)
        
        if length(F_coefficients) == 0
            break
        end
        update_F(sets_vca, F_coefficients, F_evaluations)
        
        degree += 1
    end
    
    X_train_transformed = V_to_matrix(sets_vca)
    if size(X_train_transformed, 2) != 0
        X_train_transformed = abs.(X_train_transformed)
    else
        X_train_transformed = nothing
    end
    
    return X_train_transformed, sets_vca
end


"""
Performs FindRangeNull (using SVD) for VCA.

# Arguments
- 'F::Union{Matrix{Float64}, Matrix{Int64}}': TBD
- 'C::Union{Matrix{Float64}, Matrix{Int64}}': TBD
- 'psi::Float64': vanishing parameter

# Returns
- 'V_coefficient_vectors::Union{Matrix{Float64}, Matrix{Int64}}': Coefficient vectors of polynomials we append to V.
- 'V_evaluation_vectors::Union{Matrix{Float64}, Matrix{Int64}}': Evaluation vectors of polynomials we append to V.
- 'F_coefficient_vectors::Union{Matrix{Float64}, Matrix{Int64}}': Coefficient vectors of polynomials we append to F.
- 'F_evaluation_vectors::Union{Matrix{Float64}, Matrix{Int64}}': Evaluation vectors of polynomials we append to F.
"""
function find_range_null_vca(F::Union{Matrix{Float64}, Matrix{Int64}}, C::Union{Matrix{Float64}, Matrix{Int64}}, psi::Float64)
    tmp_coefficient_vectors = zeros(Float64, size(F, 2) + size(C, 2), 0)
    tmp_evaluation_vectors = zeros(Float64, size(C, 1), 0)
    for i in 1:size(C, 2)
        tmp_coefficient_vector = zeros(1, size(C, 2))
        tmp_coefficient_vector[:, i] .= 1
        tmp_evaluation_vector = C[:, i]
        orthogonal_components = orthogonal_projection(F, tmp_evaluation_vector)
        tmp_coefficient_vector = hcat(-orthogonal_components', tmp_coefficient_vector)
        tmp_coefficient_vector = reshape(tmp_coefficient_vector, 
        size(tmp_coefficient_vector, 2), size(tmp_coefficient_vector, 1))        
        
        tmp_evaluation_vector = tmp_evaluation_vector - (F * orthogonal_components)
        
        @assert all(abs.(hcat(F, C) * tmp_coefficient_vector - tmp_evaluation_vector) .<= 1.0e-8) "sanity check"
        
        tmp_coefficient_vectors = hcat(tmp_coefficient_vectors, tmp_coefficient_vector)
        tmp_evaluation_vectors = hcat(tmp_evaluation_vectors, tmp_evaluation_vector)
    end
    
    A = tmp_evaluation_vectors
    if size(A, 1) > size(A, 2)
        A_squared = A' * A
        L, D, U = svd(A_squared; full=true)
    else
        L, D, U = svd(A; full=true)
    end
    U = U'
    
    V_coefficient_vectors = zeros(Float64, size(tmp_coefficient_vectors, 1), 0)
    V_evaluation_vectors = zeros(Float64, size(A, 1), 0)
    F_coefficient_vectors = zeros(Float64, size(tmp_coefficient_vectors, 1), 0)
    F_evaluation_vectors = zeros(Float64, size(A, 1), 0)
    
    for i in 1:size(C, 2)
        row_U = U[i, :]
        coefficient_vector = tmp_coefficient_vectors * row_U
        evaluation_vector = A * row_U
        
        loss = 1/size(A, 1) * norm(evaluation_vector, 2)^2
        
        if loss > psi
            nrm = norm(evaluation_vector, 2)
            F_coefficient_vectors = hcat(F_coefficient_vectors, (1/nrm * coefficient_vector))
            F_evaluation_vectors = hcat(F_evaluation_vectors, (1/nrm * evaluation_vector))
        else
            V_coefficient_vectors = hcat(V_coefficient_vectors, coefficient_vector)
            V_evaluation_vectors = hcat(V_evaluation_vectors, evaluation_vector)
        end
    end

    return V_coefficient_vectors, V_evaluation_vectors, F_coefficient_vectors, F_evaluation_vectors    
end


"""Evaluates transformation corresponding to the polynomials in V."""
function evaluate_vca(sets::SetsVCA, X_test::Matrix{Float64})
    X_test_transformed, sets_VCA_test = apply_V_transformation(sets, X_test)
    if size(X_test_transformed, 2) > 0
        X_test_transformed = abs.(X_test_transformed)
    else
        X_test_transformed = nothing
    end
    return X_test_transformed, sets_VCA_test
end
