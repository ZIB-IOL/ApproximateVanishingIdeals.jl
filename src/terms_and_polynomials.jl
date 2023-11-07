using LinearAlgebra

"""
Creates and keeps track of sets O and G for OAVI.
"""
mutable struct SetsOandG
    
# border terms before purging by degree
border_terms_raw
# evaluations of border_terms_raw over X by degree
border_evaluations_raw
# border terms after purging by degree
border_terms_purged
# evaluations of border_terms_purged over X by degree
border_evaluations_purged
# indices such that border_terms_raw[i][non_purging_indices[i], :] = border_terms_purged[i]
non_purging_indices   

# array of O terms by degree
O_terms
# array of evaluation of O terms over X by degree
O_evaluations
# indices in the border that get appended to O
O_indices
# indices where terms of degrees start. O_degree_indices[i] = first term of degree i
O_degree_indices

# sets for vanishing polynomials
G_coefficient_vectors
G_evaluations

# leading terms
leading_terms
    
end


"""initializes SetsOandG instance w.r.t. X_train"""
function construct_SetsOandG(X_train)
    m, n = size(X_train)
    
    # border sets
    border_terms_raw = Vector{Any}([nothing])
    border_evaluations_raw = Vector{Any}([nothing])
    border_terms_purged = Vector{Any}([nothing])
    border_evaluations_purged = Vector{Any}([nothing])
    non_purging_indices = Vector{Any}([nothing])
    
    # O sets
    O_terms = zeros(Int64, n, 1)
    O_evaluations = ones(Float64, m, 1)
    O_indices = []
    O_degree_indices = [2] # degree 1 always starts at index 2
    
    # G sets
    G_coefficient_vectors = Vector{Any}([nothing])
    G_evaluations = zeros(Float64, m, 0)
    
    # leading terms
    leading_terms = nothing
    
    return SetsOandG(border_terms_raw, border_evaluations_raw, border_terms_purged, border_evaluations_purged, non_purging_indices,
                    O_terms, O_evaluations, O_indices, O_degree_indices,
                    G_coefficient_vectors, G_evaluations,
                    leading_terms)
end


"""
updates border sets
"""
function update_border(sets, border_terms_raw, border_evaluations_raw, non_purging_indices)
    sets.non_purging_indices = append!(sets.non_purging_indices, [non_purging_indices])
    sets.border_terms_raw = append!(sets.border_terms_raw, [border_terms_raw])
    sets.border_evaluations_raw = append!(sets.border_evaluations_raw, [border_evaluations_raw])
    sets.border_terms_purged = append!(sets.border_terms_purged, [border_terms_raw[:, non_purging_indices]])
    sets.border_evaluations_purged = append!(sets.border_evaluations_purged, [border_evaluations_raw[:, non_purging_indices]])
end    


"""
updates O sets
"""
function update_O(sets, O_terms, O_evaluations, O_indices)
    sets.O_terms = hcat(sets.O_terms, O_terms)
    sets.O_evaluations = hcat(sets.O_evaluations, O_evaluations)
    sets.O_indices = append!(sets.O_indices, O_indices)
end


"""
Appends polynomial with coefficient vector based on coefficient_vector and term to G_coefficient_vectors.
"""
function update_coefficient_vectors(G_coefficient_vectors, coefficient_vector; first=false)
    if G_coefficient_vectors == nothing
        G_coefficient_vectors = coefficient_vector
    else
        if first
            lt_indices = find_first_non_zero_entries(G_coefficient_vectors)
        else
            lt_indices = find_last_non_zero_entries(G_coefficient_vectors)
        end
        
        removable_set = Set(lt_indices)
        
        indices = [x for x in 1:size(G_coefficient_vectors, 1) if x ∉ removable_set]
        
        if length(indices) == size(coefficient_vector, 1)
            updated_coefficient_vector = zeros(size(G_coefficient_vectors, 1), 1)
            updated_coefficient_vector[indices, :] = coefficient_vector
            G_coefficient_vectors = hcat(G_coefficient_vectors, updated_coefficient_vector)
        end
    end
    return G_coefficient_vectors
end


"""
updates G sets
"""
function update_G(sets, G_coefficient_vectors=nothing)
    if G_coefficient_vectors != nothing
        if size(sets.G_evaluations, 2) == 0
            sets.G_evaluations = hcat(sets.O_evaluations, sets.border_evaluations_purged[end]) * G_coefficient_vectors
        else
            current_G_evaluations = hcat(sets.O_evaluations, sets.border_evaluations_purged[end]) * G_coefficient_vectors
            sets.G_evaluations = hcat(sets.G_evaluations, current_G_evaluations)
        end
    end
    sets.G_coefficient_vectors = append!(sets.G_coefficient_vectors, [G_coefficient_vectors])
end


"""
updates leading terms
"""
function update_leading_terms(sets, leading_terms=nothing)
    if sets.leading_terms == nothing
        if leading_terms != nothing
            sets.leading_terms = leading_terms
        end
    else
        sets.leading_terms = hcat(sets.leading_terms, leading_terms)
    end
end


"""applies the transformation corresponding to G to X_test"""
function apply_G_transformation(sets::SetsOandG, X_test)
    m, n = size(X_test)
    test_sets_avi = construct_SetsOandG(X_test)
    
    append!(test_sets_avi.border_evaluations_raw, [X_test])
    
    # degree-1 border
    i = 1
    append!(test_sets_avi.border_evaluations_purged, [X_test])
    update_G(test_sets_avi, sets.G_coefficient_vectors[i+1])
    if i < size(sets.O_evaluations, 2)
        test_sets_avi.O_evaluations = hcat(test_sets_avi.O_evaluations, X_test[:, sets.O_indices[i]])
    end
    
    # higher degree border
    if sets.O_degree_indices != []  # only if higher degrees need to be considered
        i = 2
        while i < max(sets.O_degree_indices...)
            display(test_sets_avi.O_evaluations)
            display(sets.O_degree_indices)
            display(sets.O_evaluations)
            border_test_purged = reconstruct_border(test_sets_avi.border_evaluations_raw[2], 
            test_sets_avi.O_evaluations[:, sets.O_degree_indices[i]:end],
            sets.non_purging_indices[i])  

            append!(test_sets_avi.border_evaluations_purged, [border_test_purged])

            update_G(test_sets_avi, sets.G_coefficient_vectors[i+1])
            if i < size(sets.O_evaluations, 2)
                test_sets_avi.O_evaluations = hcat(test_sets_avi.O_evaluations, border_test_purged[:, sets.O_indices[i]])
            end

            i +=1
        end
    end
    return test_sets_avi.G_evaluations, test_sets_avi    
end


"""reconstructs border for O_test"""
function reconstruct_border(O1_test::Matrix{Float64}, O_test::Matrix{Float64}, non_purging_indices::Vector{Int64}; 
        X_test::Matrix{Float64}=zeros(Float64, 0, 0))
    if size(X_test) == (0, 0)
        O1_test_tile = tile(O1_test, size(O_test, 2))[:, non_purging_indices]
        O_test_repeat = repeat(O_test', outer=size(O1_test, 2))'[:, non_purging_indices]
        border_test = O1_test_tile .* O_test_repeat
    else
        border_test = X_test
    end
    return border_test
end


"""
constructs the border of 'terms'

# Arguments
- 'terms::Matrix{Int64}': Matrix with monomial terms as columns
- 'terms_evaluated::Matrix{Float64}': Matrix with evaluations of 'terms' over X
- 'X_train::Vector{Vector{Float64}}': data
- 'degree_1_terms::Matrix{Int64}': Matrix with degree 1 monomials as columns
- 'degree_1_terms_evaluated::Matrix{Float64}': evaluations of 'degree_1_terms' over X
- 'purging_terms::Matrix{Int64}': purge terms in 'terms' divisible by any of these 

# Returns
- 'border_terms_raw::Matrix{Int64}': non-purged border constructed from 'terms'
- 'border_evaluations_raw::Matrix{Float64}': non-purged evaluations of border terms over X
- 'non_purging_indices::Vector{Int64}': array of non-purging indices
"""
function construct_border(terms::Matrix{Int64}, terms_evaluated::Matrix{Float64}, X_train::Union{Matrix{Float64}, Vector{Vector{Float64}}}, 
    degree_1_terms::Matrix{Int64}=zeros(Int64, 0, 0), degree_1_terms_evaluated::Matrix{Float64}=zeros(Float64, 0, 0), 
    purging_terms::Matrix{Int64}=zeros(Int64, 0, 0))

    if typeof(X_train) != Matrix{Float64}
        X_train = vecvec_to_mat(X_train)
    end
        
    dim = length(X_train[1, :])
    
    if size(degree_1_terms, 2) == 0
        
        border_terms_raw = 1 * Matrix(I, dim, dim)
        border_evaluations_raw = X_train
    
    else
        # terms
        len_deg1 = size(degree_1_terms, 2)
        len_terms = size(terms, 2)
        
        terms_repeat = transpose(repeat(terms', outer=len_deg1))
        degree_1_tile = tile(degree_1_terms, len_terms)
        
        border_terms_raw = degree_1_tile + terms_repeat
        
        # evaluations
        len_deg1_eval = size(degree_1_terms_evaluated, 2)
        len_terms_eval = size(terms_evaluated, 2)
        
        terms_evaluated_repeat = transpose(repeat(terms_evaluated', outer=len_deg1_eval))
        degree_1_evaluated_tile = tile(degree_1_terms_evaluated, len_terms_eval)
        
        border_evaluations_raw = degree_1_evaluated_tile .* terms_evaluated_repeat
        
    end
    
    border_terms_raw, border_evaluations_raw, _ = deg_lex_sort(border_terms_raw, border_evaluations_raw)
    
    border_terms_purged, border_evaluations_purged, unique_indices = get_unique_columns(
        border_terms_raw, border_evaluations_raw)
    
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
    
    return border_terms_raw, border_evaluations_raw, non_purging_indices
    
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
