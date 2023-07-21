using LinearAlgebra
using SparseArrays

"""
Constructs the border of terms. 
    
# Arguments
- 'terms::Vector{Vector{Int}}': Array of terms to construct border from
- 'terms_evaluated::Vector{Vector{Float64}}': Array of evaluations of terms
- 'X_train::Vector{Vector{Float64}}': Array of n-dimensional data points
- 'degree_1_terms::Vector{Vector{Int64}}': Array of terms of degree 1 to construct border terms, Optional (Default is [])
- 'degree_1_terms_evaluated::Vector{Float64}': Array of evaluations of degree_1_terms, Optional (Default is [])
- 'purging_terms::Vector{Vector{Int64}}': Purge terms divisible by these terms, Optional (Default is [])
    
# Returns
- 'border_terms_raw::Vector{Vector{Int64}}': Array of unpurged border terms
- 'border_evaluations_raw::Vector{Vector{Float64}}': Array of evaluations of border terms
- 'non_purging_indices::Vector{Int64}': Array unique and non-purging indices in border_terms_raw
"""
function construct_border(
        terms::Vector{Vector{Int64}}, terms_evaluated::Vector{Vector{Float64}}, X_train::Vector{Vector{Float64}}, 
        degree_1_terms::Vector{Vector{Int64}}=[],  
        degree_1_terms_evaluated::Vector{Vector{Float64}}=[], 
        purging_terms::Vector{Vector{Int64}}=[])

    
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
function purge(
        terms::Vector{Vector{Int64}}, 
        terms_evaluated::Vector{Vector{Float64}}, 
        purging_terms::Vector{Vector{Int64}})
        
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
