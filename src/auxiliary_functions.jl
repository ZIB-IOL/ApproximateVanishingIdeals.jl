using LinearAlgebra

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
function get_unique_elements(x_1::Vector, x_2::Vector=[])

    sorted_x_1, sorted_x_2, sorted_list = deg_lex_sort(x_1, x_2)
    unique_indices = unique(i -> sorted_x_1[i], 1:length(sorted_x_1))
    x_1_unique = sorted_x_1[unique_indices]
    if x_2 != []
        x_2_unique = sorted_x_2[unique_indices]
    else
        x_2_unique = x_2
    end
        
    return x_1_unique, x_2_unique, unique_indices
    
end


"""
Transforms an object of type Matrix{Any} to an Array of Arrays where each row is an individual Array.
Mainly for testing, since rand(a:b, x, y) objects are of type Matrix. Probably not really needed but here nonetheless
    
# Arguments
- 'A::Matrix{Any}': the matrix to transform
- 'col_is_row::Int': If 1 (default 0), instead the columns of the matrix become individual Arrays.
    
# Returns
- 'transformed_A::Vector{Vector{Any}}': transformed Array
"""
function mat_to_arr_of_arrs(A::Matrix{Any}, col_is_row::Int64=0)

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


"""
Evaluates monomial m at point x.
    
# Arguments
- 'm::Vector{Int64}': monomial under which to evaluate x
- 'x::Vector{Float64}': point at which to evaluate m
- 'coeff=1': coefficient for monomial m, Optional
    
# Returns 
- 'evaluation::Float64': evaluation of m at point x
"""
function monomial_evaluation(m::Vector{Int64}, x::Vector{Float64}, coeff=1)
    result = 1
    
    for i in 1:length(m)
        result *= x[i]^m[i]
    end
    
    evaluation = coeff * result
    return evaluation
end


"""
Evaluates a set of points X under a single monomial m.
    
# Arguments
- 'm::Vector{Int64}': monomial under which to evaluate X
- 'X::Vector{Vector{Float64}}': set of points to evaluate
- 'coeff=1': coefficient for monomial m, Optional
    
# Returns
'evaluated::Vector{Float64}': Array of evaluations of m, evaluated[i] is the evaluation of m at point X[i].
"""    
function monomial_evaluation_set(m::Vector{Int64}, X::Vector{Vector{Float64}}, coeff=1)
    results = ones(length(X))
        
    for i in 1:length(X)
        results[i] = monomial_evaluation(m, X[i])
    end
    evaluated = coeff * results
    return evaluated
    
end


"""
Evaluates a set of points X under a set of monomials M.
    
# Arguments
- 'M::Vector{Vector{Int64}}': set of monomials under which to evaluate X
- 'X::Vector{Vector{Float64}}': set of points at which to evaluate M
- 'coeffs=1': coefficient(s) for all (each) monomial in M, Optional
    
# Returns
- 'all_evaluated::Vector{Vector{Float64}}': Array of Arrays where all_evaluated[i][j] is monomial M[i] evaluated at point X[j]
"""      
function monomial_set_evaluation_set(M::Vector{Vector{Int64}}, X::Vector{Vector{Float64}}, coeffs=1)  
    result = [[] for _ in 1:length(M)]
    
    for i in 1:length(M)
        result[i] = monomial_evaluation_set(M[i], X)
    end
    
    all_evaluated = coeffs .* result
    return all_evaluated
    
end


"""
Computes the sum of the row entries, which for matrices representing terms is equal to the degree.
    
# Arguments
- 'matrix::Vector{Vector{Int64}}': matrix whose row degrees have to be computed.
    
# Returns
- 'degree_array::Vector{Int64}': Array containing degrees of all rows.
"""
function compute_degree(matrix::Vector{Vector{Int64}})
    degree_array = sum.(matrix)
    return degree_array
end

    
"""
Sorts the rows of matrix_1 degree-lexicographically and matrix_2 accordingly
    
# Arguments
- 'matrix_1::Vector{Vector{Int64}}': term matrix to be sorted
- 'matrix_2': matrix getting sorted the same way matrix_1 is (Default is [])
    
# Returns
- 'matrix_1::Vector{Vector{Int64}}': matrix_1 sorted degree_lexicographically
- 'matrix_2::Vector{Vector{Int64}}': matrix_2 sorted like matrix_1
"""
function deg_lex_sort(matrix_1::Vector{Vector{Int64}}, matrix_2=[])
    degrees = compute_degree(matrix_1)
    for i in 1:length(degrees)
        # using reverse + Base.sortperm sorts first by degree and then 
        # from last to first variable exponent.
        matrix_1[i] = append!([degrees[i]], reverse(matrix_1[i]))
    end
    # find the permutation to correctly sort matrix_1
    sorted_list = sortperm(matrix_1)
    matrix_1 = matrix_1[sorted_list]
    # deletes first entry (degree) in each term
    matrix_1 = [reverse(deleteat!(matrix_1[i], 1)) for i in 1:length(matrix_1)]
    if matrix_2 != []
        matrix_2 = matrix_2[sorted_list]
    end
    return matrix_1, matrix_2, sorted_list
end
