using LinearAlgebra
using SparseArrays


function construct_border(
        terms::Vector{Vector{Int64}}, terms_evaluated::Vector{Vector{Float64}}, X_train::Vector{Vector{Float64}}, 
        degree_1_terms::Vector{Vector{Int64}}=[],  
        degree_1_terms_evaluated::Vector{Vector{Float64}}=[], 
        purging_terms::Vector{Vector{Int64}}=[])
    """
    Constructs the border of terms. 
    
    # Arguments
    - 'terms::Vector{Vector{Int}}': Array of terms to construct border from
    - 'terms_evaluated::Vector{Vector{Float64}}': Array of evaluations of terms
    - 'X_train::Vector{Vector{Float64}}': Array of n-dimensional data points
    - 'degree_1_terms::Vector{Vector{Int64}}': Array of terms of degree 1 to construct border terms, Optional 
        (Default is [])
    - 'degree_1_terms_evaluated::Vector{Float64}': Array of evaluations of degree_1_terms, Optional 
        (Default is [])
    - 'purging_terms::Vector{Vector{Int64}}': Purge terms divisible by these terms, Optional 
        (Default is [])
    
    # Returns
    - 'border_terms_raw::Vector{Vector{Int64}}': Array of unpurged border terms
    - 'border_evaluations_raw::Vector{Vector{Float64}}': Array of evaluations of border terms
    - 'non_purging_indices::Vector{Int64}': Array unique and non-purging indices in border_terms_raw
    """
    
    if degree_1_terms == []
        # 'I' comes from the LinearAlgebra package
        border_terms_raw = mat_to_arr_of_arrs(
            1 * Matrix(I, length(X_train[1]), length(X_train[1])))
            
        border_terms_raw_evaluation = X_train
        
    else
        # terms repeat + tile part
        terms_repeat = repeat(terms, inner=length(degree_1_terms))    
        l_terms = length(terms_repeat)
        m_terms = length(degree_1_terms)
        k_terms = Int(l_terms / m_terms) 
        degree_1_terms_tile = repeat(degree_1_terms, outer=k_terms)        
        border_terms_raw = degree_1_terms_tile + terms_repeat
        
        # evaluation repeat + tile part
        terms_evaluated_repeat = repeat(terms_evaluated, inner=length(degree_1_terms_evaluated))
        l_eval = length(terms_evaluated_repeat)
        m_eval = length(degree_1_terms_evaluated)
        k_eval = Int(l_eval / m_eval)
        degree_1_terms_evaluated_tile = repeat(degree_1_terms_evaluated, outer=k_eval)
        border_evaluations_raw = [degree_1_terms_evaluated_tile[i] .* terms_evaluated_repeat[i] for i in 1:l_eval]
         
    end
    
    border_terms_purged, border_evaluations_purged, unique_indices = get_unique_elements(
        border_terms_raw, border_evaluations_raw)
    
    if purging_terms != []
        border_terms_purged, border_evaluations_purged, unique_indices_2 = purge(
            border_terms_purged, border_evaluations_purged, purging_terms)
        
        if unique_indices_2 != []
            non_purging_indices = [unique_indices[i] for i in unique_indices_2]
        else
            non_purging_indices = unique_indices
        end
    else
        non_purging_indices = unique_indices
    end
    
    return border_terms_raw, border_evaluations_raw, non_purging_indices
    
end


function purge(
        terms::Vector{Vector{Int64}}, 
        terms_evaluated::Vector{Vector{Float64}}, 
        purging_terms::Vector{Vector{Int64}})
    """
    Purges all the terms in 'terms' that are divisible by at least one term in 'purging_terms'.
    
    # Arguments
    - 'terms::Vector{Vector{Int64}}': Array of terms 
    - 'terms_evaluated:::Vector{Vector{Float64}}': Array of evaluations of terms 
    - 'purging_terms::Vector{Vector{Int64}}':  Array of possible divisor terms
    
    # Returns
    - 'terms[indices]::Vector{Vector{Int64}}': purged version of terms
    - 'terms_evaluated[indices]::Vector{Vector{Float64}}': purged version of terms_evaluated
    - 'indices::Vector{Int64}': Array of indices of purged terms in 'terms'
    """
    indices = [x for x in 1:length(terms)]
    purge_indices = []
    for i in eachindex(terms)
        for j in eachindex(purging_terms)
            if all(terms[i] .- purging_terms[j] .>= 0)
                append!(purge_indices, i)
                break
            end
        end
    end
    indices = deleteat!(indices, purge_indices)
    return terms[indices], terms_evaluated[indices], indices
end


function get_unique_elements(x_1::Vector, x_2::Vector=[])
    
    """
    Finds indices of unique elements in Array x.
    
    # Arguments
    - 'x_1::Vector{Any}': object of which to find unique indices
    - 'x_2::Vector{Any}': Array with same length as x_1, Optional (Default is [])
    
    # Returns
    - 'x_1_unique::Vector{Any}': Array of unique entries in x_1
    - 'x_2_unique::Vector{Any}': Array of entries corresponding to unique elements in x_1
    - 'unique_indices::Vector{Int}': Array with positions of unique elements in x
    """
    unique_indices = unique(i -> x_1[i], 1:length(x_1))
    x_1_unique = x_1[unique_indices]
    if x_2 != []
        x_2_unique = x_2[unique_indices]
    else
        x_2_unique = x_2
    end
        
    return x_1_unique, x_2_unique, unique_indices
    
end


function mat_to_arr_of_arrs(A::Matrix{Any}, col_is_row::Int64=0)
    """
    Transforms an object of type Matrix{Any} to an Array of Arrays where each row is an individual Array.
    Mainly for testing, since rand(a:b, x, y) objects are of type Matrix. Probably not really needed but here nonetheless
    
    # Arguments
    - 'A::Matrix{Any}': the matrix to transform
    - 'col_is_row::Int': If 1 (default 0), instead the columns of the matrix become individual Arrays.
    
    # Returns
    - 'transformed_A::Vector{Vector{Any}}': transformed Array
    """

    if  col_is_row == 0
        # converts matrix to array of arrays with rows of matrix being the individual entries of array,
        # i.e. first row of matrix A becomes the array result[1]
        transformed_A = [A[i,:] for i in 1:size(A,1)]
    
    
    elseif col_is_row == 1
        # converts matrix to array of arrays with columns of matrix being the individual entries of array,
        # i.e. first column of matrix A becomes the array result[1]
        transformed_A = [A[:,i] for i in 1:size(A,2)]
       
    end
    
    return transformed_A
    
end


function monomial_evaluation(m::Vector{Int64}, x::Vector{Float64}, coeff=1)
    """
    Evaluates monomial m at point x.
    
    # Arguments
    - 'm::Vector{Int64}': monomial under which to evaluate x
    - 'x::Vector{Float64}': point at which to evaluate m
    - 'coeff=1': coefficient for monomial m, Optional
    
    # Returns 
    - 'evaluation::Float64': evaluation of m at point x
    """
    result = 1
    
    for i in 1:length(m)
        result *= x[i]^m[i]
    end
    
    evaluation = coeff * result
    return evaluation
end


function monomial_evaluation_set(m::Vector{Int64}, X::Vector{Vector{Float64}}, coeff=1)
    """
    Evaluates a set of points X under a single monomial m.
    
    # Arguments
    - 'm::Vector{Int64}': monomial under which to evaluate X
    - 'X::Vector{Vector{Float64}}': set of points to evaluate
    - 'coeff=1': coefficient for monomial m, Optional
    
    # Returns
    'evaluated::Vector{Float64}': Array of evaluations of m, evaluated[i] is the evaluation of m at point X[i].
    """    
    results = ones(length(X))
        
    for i in 1:length(X)
        results[i] = monomial_evaluation(m, X[i])
    end
    evaluated = coeff * results
    return evaluated
    
end


function monomial_set_evaluation_set(M::Vector{Vector{Int64}}, X::Vector{Vector{Float64}}, coeffs=1)
    """
    Evaluates a set of points X under a set of monomials M.
    
    # Arguments
    - 'M::Vector{Vector{Int64}}': set of monomials under which to evaluate X
    - 'X::Vector{Vector{Float64}}': set of points at which to evaluate M
    - 'coeffs=1': coefficient(s) for all (each) monomial in M, Optional
    
    # Returns
    - 'all_evaluated::Vector{Vector{Float64}}': Array of Arrays where all_evaluated[i][j] is 
    monomial M[i] evaluated at point X[j]
    """        
    result = [[] for _ in 1:length(M)]
    
    for i in 1:length(M)
        result[i] = monomial_evaluation_set(M[i], X)
    end
    
    all_evaluated = coeffs .* result
    return all_evaluated
    
end
