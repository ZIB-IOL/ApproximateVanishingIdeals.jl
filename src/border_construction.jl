"""reconstructs border for O_test"""
function reconstruct_border(O1_test::Matrix{Float64}, 
                            O_test::Matrix{Float64}, 
                            non_purging_indices::Vector{Int64},
                            permutation::Vector{Int64}; 
                            X_test::Matrix{Float64}=zeros(Float64, 0, 0)
                            )
    if size(X_test) == (0, 0)
        O1_test_tile = repeat(O1_test, inner=[1, size(O_test, 2)])
        O_test_repeat = repeat(O_test, outer=[1, size(O1_test, 2)])
        border_test = O1_test_tile .* O_test_repeat
        border_test = border_test[:, permutation]
    else
        border_test = X_test
    end
    return border_test[:, non_purging_indices]
end


"""
constructs the border of 'terms'

# Arguments
- 'terms::Matrix{Int64}': Matrix with monomial terms as columns
- 'terms_evaluated::Matrix{Float64}': Matrix with evaluations of 'terms' over X
- 'X_train::Matrix{Float64}': data
- 'degree_1_terms::Matrix{Int64}': Matrix with degree 1 monomials as columns
- 'degree_1_terms_evaluated::Matrix{Float64}': evaluations of 'degree_1_terms' over X
- 'purging_terms::Matrix{Int64}': purge terms in 'terms' divisible by any of these 

# Returns
- 'border_terms_raw::Matrix{Int64}': non-purged border constructed from 'terms'
- 'border_evaluations_raw::Matrix{Float64}': non-purged evaluations of border terms over X
- 'non_purging_indices::Vector{Int64}': array of non-purging indices
- 'raw_permutation::Vector{Int64}': array with deg-lex ordering permutation
"""
function construct_border(terms::Matrix{Int64}, 
                          terms_evaluated::Matrix{Float64}, 
                          X_train::Matrix{Float64}, 
                          degree_1_terms::Matrix{Int64}=zeros(Int64, 0, 0), 
                          degree_1_terms_evaluated::Matrix{Float64}=zeros(Float64, 0, 0), 
                          purging_terms::Matrix{Int64}=zeros(Int64, 0, 0)
                          )
        
    dim = size(X_train, 2)  
    
    if size(degree_1_terms, 2) == 0
        
        border_terms_raw = 1 * Matrix(I, dim, dim)
        border_evaluations_raw = X_train
    
    else        
        # tile and repeat for creating next degree monomials
        terms_repeat = repeat(terms, outer=[1, size(degree_1_terms, 2)])
        degree_1_tile = repeat(degree_1_terms, inner=[1, size(terms, 2)])
    
        border_terms_raw = degree_1_tile .+ terms_repeat
        
        # tile and repeat for creating next degree monomial evaluations
        terms_evaluated_repeat = repeat(terms_evaluated, outer=[1, size(degree_1_terms_evaluated, 2)])
        degree_1_evaluated_tile = repeat(degree_1_terms_evaluated, inner=[1, size(terms_evaluated, 2)])
        
        border_evaluations_raw = degree_1_evaluated_tile .* terms_evaluated_repeat
        
    end
    
    border_terms_raw, border_evaluations_raw, raw_permutation = deg_lex_sort(border_terms_raw, border_evaluations_raw)
    
    border_terms_purged, border_evaluations_purged, unique_indices = get_unique_columns(border_terms_raw, border_evaluations_raw)
    
    if size(purging_terms, 2) != 0
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
    
    return border_terms_raw, border_evaluations_raw, non_purging_indices, raw_permutation
    
end


"""
purges each term in 'terms' that is divisible by at least one term 'purging_terms'

# Arguments
- 'terms::Matrix{Int64}': Matrix with monomial terms as columns
- 'terms_evaluated::Matrix{Float64}': evaluations of 'terms' over data
- 'purging_terms::Matrix{Int64}': Matrix with purging terms as columns

# Returns
- 'terms[:, inidces]::Matrix{Int64}': purged version of terms
- 'terms_evaluated[:, indices]::Matrix{Float64}': purged evaluations
- 'indices::Vector{Int64}': array with non-purging indices
"""
function purge(terms::Matrix{Int64}, terms_evaluated::Matrix{Float64}, purging_terms::Matrix{Int64})
    purging_indices = []
    indices = [x for x in 1:size(terms, 2)]
    
    for i in 1:size(terms, 2)
        for j in 1:size(purging_terms, 2)
            if all(terms[:, i] - purging_terms[:, j] .>= 0)
                append!(purging_indices, i)
                break
            end
        end
    end
    
    indices = deleteat!(indices, purging_indices)
    return terms[:, indices], terms_evaluated[:, indices], indices
end
