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
    append!(sets.O_indices, [O_indices])
end


"""
updates G sets
"""
function update_G(sets, G_coefficient_vectors=nothing)
    if G_coefficient_vectors !== nothing
        if size(sets.G_evaluations, 2) === 0
            sets.G_evaluations = hcat(sets.O_evaluations, sets.border_evaluations_purged[end]) * G_coefficient_vectors
        else
            current_G_evaluations = hcat(sets.O_evaluations, sets.border_evaluations_purged[end]) * G_coefficient_vectors
            sets.G_evaluations = hcat(sets.G_evaluations, current_G_evaluations)
        end
    end
    append!(sets.G_coefficient_vectors, [G_coefficient_vectors])
end


"""
updates leading terms
"""
function update_leading_terms(sets, leading_terms=nothing)
    if sets.leading_terms === nothing
        if leading_terms !== nothing
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
    if sets.O_degree_indices !== []  # only if higher degrees need to be considered
        i = 2
        while i < max(sets.O_degree_indices...)
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


"""
Evaluates the transformation corresponding to the polynomials in G_coefficient_vectors.

# Arguments
- 'sets::SetsOandG': instance of SetsOandG, containing information on transformation

# Returns
- 'total_number_of_zeros::Int64': Sum of all zero entries in coefficient_vectors in G_coefficient_vectors.
- 'total_number_of_entries::Int64': Total number of entries in coefficient_vectors in G_coefficient_vectors.
- 'avg_sparsity::Float64': The average sparsity of coefficient_vectors in G_coefficient_vectors.
- 'number_of_polynomials::Int64': Number of polynomials in G.
- 'number_of_terms::Int64': Number of terms in O.
- 'degree::Float64': Average degree of polynomials in G.
"""
function evaluate_transformation(sets::SetsOandG)
    number_of_polynomials, total_number_of_zeros, total_number_of_entries, avg_sparsity, degree = 0, 0, 0, 0.0, 0
    for i in 1:size(sets.G_coefficient_vectors, 1)
        coefficient_vectors = sets.G_coefficient_vectors[i]
        if coefficient_vectors != nothing && size(coefficient_vectors, 2) > 0
            degree += (i-1) * size(coefficient_vectors, 2)  # (i-1) because we start with degree 0 but julia indexing starts at 1
            indices = find_last_non_zero_entries(coefficient_vectors)
            for j in 1:size(coefficient_vectors, 2)
                poly = coefficient_vectors[1:indices[j], j]
                number_of_entries = size(poly, 1) - (j-1)  
                number_non_zeros = sum(poly .!= 0)
                number_of_zeros = max(number_of_entries - number_non_zeros, 0)
                sparsity = number_of_zeros / max(number_of_entries - 1, 1)
                total_number_of_entries += number_of_entries
                total_number_of_zeros += number_of_zeros
                avg_sparsity += sparsity
                number_of_polynomials += 1                
            end
        end
    end
    if number_of_polynomials > 0
        avg_sparsity /= number_of_polynomials
        degree /= number_of_polynomials
    end
    number_of_terms = size(sets.O_evaluations, 2)
    return total_number_of_zeros, total_number_of_entries, avg_sparsity, number_of_polynomials, number_of_terms, degree    
end
