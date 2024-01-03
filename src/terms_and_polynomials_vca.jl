using LinearAlgebra

"""
Manages the sets V, C and F for VCA.
"""
mutable struct SetsVCA
    # C
    X
    Cs
    
    # F
    Fs
    F_coefficient_vectors
    
    # V
    Vs
    V_coefficient_vectors
end


"""
Given the data X, this function constructs the intial state of 'SetsVCA'.
"""
function construct_SetsVCA(X::Union{Matrix{Float64}, Matrix{Int64}})
    Cs = []
    Fs = 1/ sqrt(size(X, 1)) * [ones(Float64, size(X, 1), 1)]
    F_coefficient_vectors = [ones(Float64, 1, 1)]
   
    Vs = []
    V_coefficient_vectors = []
    return SetsVCA(X, Cs, Fs, F_coefficient_vectors, Vs, V_coefficient_vectors)
end;


"""updates F sets"""
function update_F(sets::SetsVCA, F_coefficient_vectors::Union{Matrix{Float64}, Matrix{Int64}}, F_evaluation_vectors::Matrix{Float64})
    sets.F_coefficient_vectors = append!(sets.F_coefficient_vectors, [F_coefficient_vectors])
    sets.Fs = append!(sets.Fs, [F_evaluation_vectors])
end


"""Transforms SetsVCA.Fs into single matrix."""
function F_to_matrix(sets::SetsVCA)
    F_array = zeros(size(sets.Fs[1], 1), 0)
    for F in sets.Fs
        if size(F, 2) > 0
            F_array = hcat(F_array, F)
        end
    end
    return F_array
end


"""updates V sets."""
function update_V(sets::SetsVCA, V_coefficient_vectors::Union{Matrix{Float64}, Matrix{Int64}}, V_evaluation_vectors::Matrix{Float64})
    V_coefficients, V_evaluations = zeros(Float64, 0, 0), zeros(Float64, 0, 0)
    if size(V_coefficient_vectors, 2) > 0 && size(V_evaluation_vectors, 2) > 0
        V_coefficients = V_coefficient_vectors
        V_evaluations = V_evaluation_vectors
    end
    sets.Vs = append!(sets.Vs, [V_evaluations])
    sets.V_coefficient_vectors = append!(sets.V_coefficient_vectors, [V_coefficients])
end


"""Transforms SetsVCA.Vs into single matrix."""
function V_to_matrix(sets::SetsVCA)
    V_array = zeros(size(sets.Vs[1], 1), 0)
    for V in sets.Vs
        if size(V, 2) > 0
            V_array = hcat(V_array, V)
        end
    end
    return V_array
end


"""
Applies transformation corresponding to V polynomials to X_test.
"""
function apply_V_transformation(sets::SetsVCA, X_test::Matrix{Float64})
    sets_VCA_test = construct_SetsVCA(X_test)
    degree = 1
    while degree <= size(sets.Cs, 1)
        border = construct_border(sets_VCA_test; degree=degree)
        update_C(sets_VCA_test, border)
        degree += 1
        
        V_coefficients_test = nothing
        if degree - 1 <= size(sets.V_coefficient_vectors, 1)
            V_coefficients_test = sets.V_coefficient_vectors[degree - 1]
        end
        
        if V_coefficients_test isa Matrix || V_coefficients_test isa Vector
            stack_Fs_and_Cs = hcat(F_to_matrix(sets_VCA_test), sets_VCA_test.Cs[end])
            V_evaluations_test = (V_coefficients_test' * stack_Fs_and_Cs')'
            
            # change V_evaluations_test from 'Adjoint' type to 'Matrix' type
            V_evaluations_test = Matrix(V_evaluations_test)
        else
            V_evaluations_test = nothing
        end

        update_V(sets_VCA_test, V_coefficients_test, V_evaluations_test)
        
        F_coefficients_test = nothing
        if degree <= size(sets.F_coefficient_vectors, 1)
            F_coefficients_test = sets.F_coefficient_vectors[degree]
        end
        
        if F_coefficients_test isa Matrix ||F_coefficients_test isa Vector
            stack_Fs_and_Cs = stack_Fs_and_Cs = hcat(F_to_matrix(sets_VCA_test), sets_VCA_test.Cs[end])
            F_evaluations_test = (F_coefficients_test' * stack_Fs_and_Cs')'
            
            # change F_evaluations_test from 'Adjoint' type to 'Matrix' type
            F_evaluations_test = Matrix(F_evaluations_test)
            update_F(sets_VCA_test, F_coefficients_test, F_evaluations_test)
        else
            break
        end            
    end
    X_test_transformed = V_to_matrix(sets_VCA_test)
    return X_test_transformed, sets_VCA_test
end;


"""updates C"""
function update_C(sets::SetsVCA, vectors::Matrix{Float64})
    sets.Cs = append!(sets.Cs, [vectors])
end;


"""constructs the border for current state of the algorithm."""
function construct_border(sets::SetsVCA; degree::Int64=1)
    if degree != 1
        F1 = sets.Fs[2]
        F_current = sets.Fs[end]
        F1_tile = tile(F1, size(F_current, 2))
        F_current_repeat = repeat(F_current, outer=(1, size(F1, 2)))
        border = F1_tile .* F_current_repeat
    else
        border = sets.X
    end
    return border
end;


"""
Evaluates the transformation corresponding to the polynomials in V w.r.t. the functions in F and C

# Arguments
- 'sets::SetsVCA': instance of SetsVCA.

# Returns
- 'total_number_of_zeros::Int64': Sum of all zero entries in coefficient_vectors in V_coefficient_vectors.
- 'total_number_of_entries::Int64': Total number of entries in coefficient_vectors in V_coefficient_vectors.
- 'avg_sparsity::Float64': The average sparsity of coefficient_vectors in V_coefficient_vectors.
- 'number_of_polynomials::Int64': Number of polynomials in V_coefficient_vectors.
- 'number_of_terms::Int64': Number of non-vanishing terms.
- 'degree::Float64': Average degree of polynomials in V.
"""
function evaluate_transformation(sets::SetsVCA)
    total_number_of_zeros, total_number_of_entries, avg_sparsity, degree = 0, 0, 0.0, 0
    number_of_polynomials = 0
    for i in 1:size(sets.V_coefficient_vectors, 1)
        coefficient_vectors = sets.V_coefficient_vectors[i]
        if size(coefficient_vectors, 2) > 0
            degree += i * size(coefficient_vectors, 2)
            
            for j in 1:size(coefficient_vectors, 2)
                number_of_polynomials += 1
                poly = coefficient_vectors[:, j]
                number_of_entries = length(poly)
                number_of_zeros = number_of_entries - sum((poly .!= 0))
                sparsity = number_of_zeros / max(number_of_entries - 1, 1)
                total_number_of_zeros += number_of_zeros
                total_number_of_entries += number_of_entries
                avg_sparsity += sparsity                
            end
            
        end
    end
    avg_sparsity = avg_sparsity / number_of_polynomials
    degree = degree / number_of_polynomials
    number_of_terms = size(F_to_matrix(sets), 2)
    return (total_number_of_zeros, total_number_of_entries, avg_sparsity, number_of_polynomials, number_of_terms, degree)
end;

