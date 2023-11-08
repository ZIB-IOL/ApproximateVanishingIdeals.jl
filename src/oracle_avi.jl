using LinearAlgebra

"""
Creates OAVI feature transformation fitted to X_train

# Arguments
- 'X_train::Union{Matrix{Float64}, Vector{Vector{Float64}}}': training data
- 'max_degree::Int64': max degree of polynomials computed (default 10)
- 'psi::Float64': vanishing extent (default 0.1)
- 'epsilon::Float64': accuracy for convex optimizer (default 0.001)
- 'tau::Union{Float64, Int64}': upper bound on norm of coefficient vector
- 'lambda::Float64': regularization parameter
- 'oracle::Union{String, <:Function}': string denoting which predefined constructor to use OR constructor function.
                                        (external constructor function MUST have 'data' and 'labels' as varargs + kwargs)
- 'orcl_kwargs::Vector': Array containing keyword arguments for external constructor functions
                            e.g.: orcl_kwargs = [(:psi, 0.05), (:epsilon, 0.0003)] inputs psi=0.05 and epsilon=0.0003 as kwargs

# Returns
- 'X_train_transformed::Vector{Vector{Float64}}': transformed X_train
- 'sets::SetsOandG': instance of mutable struct 'SetsOandG' keeping track of important sets 
""" 
function fit(X_train::Union{Matrix{Float64}, Vector{Vector{Float64}}}; 
        max_degree::Int64=10, psi::Float64=0.1, epsilon::Float64=0.001, tau::Union{Float64, Int64}=1000.,
        lambda::Float64=0., oracle::Union{String, <:Function}="CG", max_iters::Int64=10000, 
        inverse_hessian_boost::String="false", orcl_kwargs::Vector=[]) 

    if typeof(X_train) != Matrix{Float64}
        X_train = vecvec_to_mat(X_train)
    end
    
    m, n = size(X_train)
    
    sets = construct_SetsOandG(X_train)
    
    degree = 0
    while degree < max_degree
        degree += 1        
        
        if degree == 1
            border_terms_raw, border_evaluations_raw, non_purging_indices = construct_border(sets.O_terms, sets.O_evaluations, X_train)
        else
            deg_idx = sets.O_degree_indices[degree]
            border_terms_raw, border_evaluations_raw, non_purging_indices = construct_border(sets.O_terms[:, deg_idx:end], sets.O_evaluations[:, deg_idx:end], X_train, 
            sets.border_terms_raw[2], sets.border_evaluations_raw[2])
        end
        
        border_terms = border_terms_raw[:, non_purging_indices]
        border_evaluations = border_evaluations_raw[:, non_purging_indices]
        
        update_border(sets, border_terms_raw, border_evaluations_raw, non_purging_indices)
        
        O_indices = []
        leading_terms = []
        G_coefficient_vectors = nothing

        data = sets.O_evaluations
        data_squared = data' * data
        data_squared_inverse = nothing
        
        if inverse_hessian_boost in ["weak", "full"]
            data_squared_inverse = inv(data_squared)
        end

        for col_idx in 1:size(border_terms, 2)

            if G_coefficient_vectors != nothing
                G_coefficient_vectors = vcat(G_coefficient_vectors, zeros(Float64, 1, size(G_coefficient_vectors, 2)))
            end

            term_evaluated = border_evaluations[:, col_idx] 
            data_term_evaluated = data' * term_evaluated
            term_evaluated_squared = term_evaluated' * term_evaluated
            
            if oracle in ["CG", "BCG", "BPCG"]
                coefficient_vector, loss = conditional_gradients(oracle, data, term_evaluated, 
                lambda, data_squared, data_term_evaluated, term_evaluated_squared; data_squared_inverse=data_squared_inverse, 
                psi=psi, epsilon=epsilon, tau=1. * tau, inverse_hessian_boost=inverse_hessian_boost)
            
            elseif oracle == "ABM"
                coefficient_vector, loss = abm(data, labels, data_squared, data_term_evaluated, term_evaluated_squared)
            
            # external constructor
            else
                coefficent_vector, loss = oracle(data, labels; orcl_kwargs...)
            end
            
            if loss <= psi
                leading_terms = append!(leading_terms, col_idx)
                G_coefficient_vectors = update_coefficient_vectors(G_coefficient_vectors, coefficient_vector)
            else
                O_indices = append!(O_indices, col_idx)
                data, data_squared, data_squared_inverse = streaming_matrix_updates(data, data_squared, data_term_evaluated,
                                                            term_evaluated, term_evaluated_squared; A_squared_inv=data_squared_inverse)
            end
            
        end
        
        update_leading_terms(sets, border_terms[:, leading_terms])
        update_G(sets, G_coefficient_vectors)
        
        if O_indices == []
            sets.O_degree_indices = sets.O_degree_indices[1:length(sets.O_degree_indices)-1]
            break
        else
            update_O(sets, border_terms[:, O_indices], border_evaluations[:, O_indices], O_indices)
            sets.O_degree_indices = append!(sets.O_degree_indices, size(sets.O_terms, 2) + 1)
        end
    end
    
    X_train_transformed = sets.G_evaluations
    
    if size(X_train_transformed, 2) != 0
        X_train_transformed = abs.(X_train_transformed)
    else
        X_train_transformed = nothing
    end
    
    return X_train_transformed, sets
    
end
