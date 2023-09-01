"""
Creates OAVI feature transformation fitted to X_train

# Arguments
- 'X_train::Vector{Vector{Float64}}': training data
- 'max_degree::Int64': max degree of polynomials computed (default 10)
- 'psi::Float64': vanishing extent (default 0.1)
- 'epsilon::Float64': accuracy for convex optimizer (default 0.001)
- 'tau::Union{Float64, Int64}': upper bound on norm of coefficient vector

# Returns
- 'X_train_transformed::Vector{Vector{Float64}}': transformed X_train
- 'sets::SetsOandG': instance of mutable struct 'SetsOandG' keeping track of important sets 
""" 
function fit(X_train::Union{Matrix{Float64}, Vector{Vector{Float64}}}; 
        max_degree::Int64=10, psi::Float64=0.1, epsilon::Float64=0.001, tau::Union{Float64, Int64}=1000,
        lmbda::Float64=0., tol::Float64=0.0001, objective_type::String="L2Loss", region_type::String="L1Ball", 
        oracle_type::String="CG", max_iters::Int64=10000, inverse_hessian_boost::String="false")

    if typeof(X_train) != Matrix{Float64}
        X_train = vecvec_to_mat(X_train)
    end
    
    m, n = size(X_train)
    
    sets = SetsOandG(Vector{Any}([nothing]), Vector{Any}([nothing]), Vector{Any}([nothing]), Vector{Any}([nothing]), Vector{Any}([nothing]),
    zeros(Int64, n, 1), ones(Float64, m, 1), [], Vector{Int64}([]),
    Vector{Any}([nothing]), zeros(Float64, m, 0), 
    nothing)
    
    degree = 0
    while degree < max_degree
        degree += 1 
        sets.O_degree_indices = append!(sets.O_degree_indices, size(sets.O_terms, 2) + 1)       
        
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
            data_with_labels = hcat(data, term_evaluated)
            
            f, grad!, region = nothing, nothing, nothing           
            
            if objective_type == "L2Loss"
                objective_data, f, grad! = construct_L2Loss(data, term_evaluated; lmbda=lmbda, data_squared=data_squared, labels_squared=term_evaluated_squared, 
                                    data_squared_inverse=data_squared_inverse, data_labels=data_term_evaluated)
            end

            if region_type == "L1Ball"
                region = FrankWolfe.LpNormLMO{1}(tau-1)
            end

            @assert f != nothing "Objective function f not defined."
            @assert grad! != nothing "Gradient of f not defined."
            @assert region != nothing "Feasible region not defined."

            if inverse_hessian_boost == "full"
                x0 = l1_projection(objective_data.solution; radius=tau-1)
                x0 = reshape(x0, size(x0, 1))
                
                coefficient_vector = call_oracle(f, grad!, region, x0)
                coefficient_vector = vcat(coefficient_vector, [1])
                
                loss = (1 / size(data, 1)) * norm(data_with_labels * coefficient_vector, 2)^2
                
            elseif inverse_hessian_boost == "weak"
                x0 = l1_projection(objective_data.solution; radius=tau-1)
                x0 = reshape(x0, size(x0, 1))
                
                coefficient_vector = call_oracle(f, grad!, region, x0; oracle=oracle_type)
                coefficient_vector = vcat(coefficient_vector, [1])
                
                loss = (1 / size(data, 1)) * norm(data_with_labels * coefficient_vector, 2)^2
                
                if loss <= psi
                    x0 = compute_extreme_point(region, zeros(Float64, size(data, 2)))
                    x0 = Vector(x0)
                    
                    tmp_coefficient_vector = call_oracle(f, grad!, region, x0; oracle=oracle_type)
                    tmp_coefficient_vector = vcat(tmp_coefficient_vector, [1])
                    
                    loss_2 = (1 / size(data, 1)) * norm(data_with_labels * tmp_coefficient_vector, 2)^2
                    
                    if loss_2 <= psi
                        loss = loss_2
                        coefficient_vector = tmp_coefficient_vector
                    end
                    
                end     
                
            else
                x0 = compute_extreme_point(region, zeros(Float64, size(data, 2)))
                x0 = Vector(x0)
                
                coefficient_vector = call_oracle(f, grad!, region, x0; oracle=oracle_type)
                coefficient_vector = vcat(coefficient_vector, [1])

                loss = (1 / size(data, 1)) * norm(data_with_labels * coefficient_vector, 2)^2
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
            break
        else
            update_O(sets, border_terms[:, O_indices], border_evaluations[:, O_indices], O_indices)
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
    
    
"""
Projects x onto the L1 Ball with radius 'radius'.

Reference: 
"Efficient Projections onto the ℓ1-Ball for Learning in High Dimensions", https://stanford.edu/~jduchi/projects/DuchiShSiCh08.pdf
"""
function l1_projection(x; radius=1.)
    @assert radius > 0 "Radius must be positive."
    
    if norm(x, 1) <= radius
        return x
    end
    
    n = size(x, 1)
    x = reshape(x, n, 1)
    v = abs.(x)
    u = sort(v, dims=1)
    u = reverse(u)
    
    csum = cumsum(u, dims=1)
    
    p = findlast(u .* collect(1:n) .> csum .- radius)[1]
    
    theta = (1 / p) * (csum[p] - radius)
    
    w = reshape([max(v[i]-theta, 0) for i in 1:size(v, 1)], n, 1)
    return sign.(x) .* w
end