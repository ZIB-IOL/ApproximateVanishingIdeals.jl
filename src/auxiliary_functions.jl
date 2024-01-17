"""
finds unique columns in matrix x1 and returns only unique elements in x1 as well as corresponding columns in x2.
"""
function get_unique_columns(x1::Matrix, x2::Matrix=zeros(0, 0))
    cx1, cx2 = copy(x1), copy(x2)
    sorted_x1, sorted_x2, sorted_list = deg_lex_sort(cx1, cx2)
    unique_indices = unique(i -> sorted_x1[:, i], 1:size(sorted_x1, 2))
    x1_unique = sorted_x1[:, unique_indices]
    if size(x2) != (0, 0)
        x2_unique = sorted_x2[:, unique_indices]
    else
        x2_unique = cx2
    end

    return x1_unique, x2_unique, unique_indices
end


"""
sorts matrix1 degree-lexicographically and matrix2 accordingly
"""
function deg_lex_sort(matrix1::Matrix, matrix2::Matrix=zeros(0, 0))
    sorted_matrix2 = nothing
    
    rev_mat1 = matrix1[end:-1:1, :]
    mat = vcat(compute_degree(matrix1), rev_mat1)
    sorted_list = sortperm(collect(eachcol(mat)))

    sorted_matrix1 = matrix1[:, sorted_list]

    if size(matrix2) != (0, 0)
        sorted_matrix2 = matrix2[:, sorted_list]
    end
    return sorted_matrix1, sorted_matrix2, sorted_list
end


"""
Computes sum (degree) of columns (terms) in A and returns output as (1 x size(A, 2)) matrix.
"""
function compute_degree(A::Matrix{Int64})
    return sum(A, dims=1)
end                     


"""
Finds last non-zero entry in each column and returns a list of the indices.
In case of a zero-column returns last index.
"""
function find_last_non_zero_entries(matrix::Matrix)
    indices = Array{Int64}([])
    @inbounds for a in eachcol(matrix)
        idx = []
        @inbounds for i in eachindex(a)
            if a[i] != 0
                push!(idx, i)
            end
        end
        if length(idx) > 0
            indices = append!(indices, idx[end])
        else
            indices = append!(indices, length(a))
        end
    end
    return indices
end
    
    
"""
Finds first non-zero entry in each column and returns a list of the indices.
In case of a zero-column returns first index.
"""
function find_first_non_zero_entries(matrix::Matrix)
    indices = Array{Int64}([])
    @inbounds for a in eachcol(matrix)
        indices = append!(indices, findmax(a)[2])
    end
    return indices
end
    

"""Obtains orthogonal projections of vector projected onto vectors."""
function orthogonal_projection(vectors::Matrix{Float64}, vector::Vector{Float64})
    @assert size(vectors, 1) == size(vector, 1) "Vectors should have identical length."
    orthogonal_components = vector' * vectors
    return orthogonal_components'
end