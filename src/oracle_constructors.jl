"""
Returns coefficient_vector and loss found through CG-based algorithm fit to 'data'.

# Arguments
- 'oracle_type::String': type of CG-based algorithm (choice from 'CG', 'BCG', 'BPCG')
- 'data::Union{Matrix{Float64}, Matrix{Int64}}': evaluations of O_terms 
- 'labels::Union{Vector{Float64}, Vector{Int64}}': current border term evaluated
- 'lambda::Union{Float64, Int64}': regularization parameter
- 'data_squared::Union{Matrix{Float64}, Matrix{Int64}}': data' * data 
- 'data_labels::Vector{Float64}': data' * labels
- 'labels_squared::Float64': labels' * labels
- 'data_squared_inverse::Union{Matrix{Float64}, Matrix{Int64}, Nothing}': inverse of data_squared for IHB, optional (default is nothing)
- 'psi::Float64': vanishing parameter, optional (default is 0.1)
- 'epsilon::Float64': solver accuracy, optional (default is 0.001)
- 'tau::Float64': bound on coefficient_vector norm, optional (default is 1000.)
- 'inverse_hessian_boost::String': whether or not to use IHB (choice from 'false', 'weak', 'full'), optional (default is 'false')

# Returns
- 'coefficient_vector::Vector{Float64}': coefficient_vector minimizing L2Loss over L1Ball
- 'loss::Float64': loss w.r.t. 'coefficient_vector'
"""
function conditional_gradients(oracle_type::String, 
        data::Union{Matrix{Float64}, Matrix{Int64}}, 
        labels::Union{Vector{Float64}, Vector{Int64}},
        lambda::Union{Float64, Int64}, 
        data_squared::Union{Matrix{Float64}, Matrix{Int64}}, 
        data_labels::Vector{Float64},
        labels_squared::Float64;
        data_squared_inverse::Union{Matrix{Float64}, Matrix{Int64}, Nothing}=nothing,
        psi::Float64=0.1,
        epsilon::Float64=0.001,
        tau::Float64=1000.,
        inverse_hessian_boost::String="false")
    
    m, n = size(data)
    data_with_labels = hcat(data, labels)
           
    solution, f, grad! = L2Loss(data, labels, lambda, data_squared, data_labels, labels_squared; data_squared_inverse=data_squared_inverse)
        
    if oracle_type == "CG"
        oracle = frank_wolfe
    elseif oracle_type == "BCG"
        oracle = blended_conditional_gradient
    elseif oracle_type == "BPCG"
        oracle = FrankWolfe.blended_pairwise_conditional_gradient
    end
    
    region = FrankWolfe.LpNormLMO{1}(tau-1)
    
    if inverse_hessian_boost in ["weak", "full"]
        x0 = l1_projection(solution; radius=tau-1)
    else
        x0 = compute_extreme_point(region, zeros(Float64, n))
        x0 = Vector(x0)
    end
    
    if inverse_hessian_boost == "weak"
        coefficient_vector, _ = oracle(f, grad!, region, x0; epsilon=epsilon)
        coefficient_vector = vcat(coefficient_vector, [1])
    
        loss = 1/m * norm(data_with_labels * coefficient_vector, 2)^2
        
        if loss <= psi
            x0 = compute_extreme_point(region, zeros(Float64, n))
            x0 = Vector(x0)
            tmp_coefficient_vector, _ = oracle(f, grad!, region, x0; epsilon=epsilon)
            tmp_coefficient_vector = vcat(tmp_coefficient_vector, [1])
            
            loss2 = 1/m * norm(data_with_labels * tmp_coefficient_vector, 2)^2
            
            if loss2 <= psi
                loss = loss2
                coefficient_vector = tmp_coefficient_vector
            end
                   
        end
    else    
        coefficient_vector, _ = oracle(f, grad!, region, x0; epsilon=epsilon)
        coefficient_vector = vcat(coefficient_vector, [1])
       
        loss = 1/m * norm(data_with_labels * coefficient_vector, 2)^2
    end        
    return coefficient_vector, loss
end


"""
Runs ABM algorithm to find coefficient vector and computes loss.

# Arguments
- 'oracle_type::String': string denoting which oracle to construct
- 'data::Union{Matrix{Float64}, Matrix{Int64}}': data (O_evaluations)
- 'labels::Union{Vector{Float64}, Vector{Int64}}': labels (term_evaluated)
- 'lambda::Union{Float64, Int64}': regularization parameter (if applicable)
- 'data_squared::Union{Matrix{Float64}, Matrix{Int64}}': squared data 
- 'data_labels::Vector{Float64}': data' * labels 
- 'labels_squared::Float64': labels' * labels
- 'data_squared_inverse::Union{Matrix{Float64}, Matrix{Int64}, Nothing}': inverse of data_squared (default is nothing)

# Returns
- 'coefficient_vector::Vector{Float64}': coefficient vector minimizing ABM optimization problem
- 'loss::Float64': loss w.r.t. 'coefficient_vector' 
"""
function abm(data::Union{Matrix{Float64}, Matrix{Int64}}, 
        labels::Union{Vector{Float64}, Vector{Int64}},
        data_squared::Union{Matrix{Float64}, Matrix{Int64}}, 
        data_labels::Vector{Float64},
        labels_squared::Float64)
    data_with_labels = hcat(data, labels)
    m = size(data_with_labels, 1)
    
    if size(data_with_labels, 1) > size(data_with_labels, 2)
        data_squared_with_labels = hcat(data_squared, data_labels)
        bottom_row = vcat(data_labels, labels_squared)
        bottom_row = bottom_row'
        data_squared_with_labels = vcat(data_squared_with_labels, bottom_row)
        F = svd(data_squared_with_labels)
    else
        F = svd(data_with_labels)
    end
    
    U, S, Vt = F.U, F.S, F.Vt
    coefficient_vector = Vt[:, end]
    loss = 1/size(data, 1) * norm(data_with_labels * coefficient_vector, 2)^2
    
    return coefficient_vector, loss
end