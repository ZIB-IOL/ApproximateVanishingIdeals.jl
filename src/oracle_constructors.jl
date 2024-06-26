"""
Returns coefficient_vector and loss found through CG-based algorithm fit to 'data'.

# Arguments
- 'oracle_type::String': type of CG-based algorithm (choice from 'CG', 'BCG', 'BPCG')
- 'data::Matrix{Float64}': evaluations of O_terms 
- 'labels::Vector{Float64}': current border term evaluated
- 'lambda::Union{Float64, Int64}': regularization parameter
- 'data_squared::Matrix{Float64}': data' * data 
- 'data_labels::Vector{Float64}': data' * labels
- 'labels_squared::Float64': labels' * labels
- 'data_squared_inverse::Union{Matrix{Float64}, Nothing}': inverse of data_squared for IHB, optional (default is nothing)
- 'psi::Float64': vanishing parameter, optional (default is 0.1)
- 'epsilon::Float64': solver accuracy, optional (default is 0.001)
- 'tau::Float64': bound on coefficient_vector norm, optional (default is 1000.)
- 'inverse_hessian_boost::String': whether or not to use IHB (choice from 'false', 'weak', 'full'), optional (default is 'false')

# Returns
- 'coefficient_vector::Vector{Float64}': coefficient_vector minimizing L2Loss over L1Ball
- 'loss::Float64': loss w.r.t. 'coefficient_vector'
"""
function conditional_gradients(
        oracle_type::String, 
        data::Matrix{Float64}, 
        labels::Vector{Float64},
        lambda::Union{Float64, Int64}, 
        data_squared::Matrix{Float64}, 
        data_labels::Vector{Float64},
        labels_squared::Float64;
        data_squared_inverse::Union{Matrix{Float64}, Nothing}=nothing,
        psi::Float64=0.1,
        epsilon::Float64=0.001,
        tau::Float64=1000.,
        inverse_hessian_boost::String="false",
        max_iters=10000)
    
    m, n = size(data)
    data_with_labels = hcat(data, labels)
    
    # Create objective function and gradient
    solution, f, grad! = L2Loss(data, labels, lambda, data_squared, data_labels, labels_squared; data_squared_inverse=data_squared_inverse)
    
    # determine oracle
    if oracle_type == "CG"
        oracle = FrankWolfe.frank_wolfe
    elseif oracle_type == "Away" || oracle_type == "AFW"
        oracle = FrankWolfe.away_frank_wolfe
    elseif oracle_type == "PCG"
        oracle = FrankWolfe.pairwise_frank_wolfe
    elseif oracle_type == "Lazy" || oracle_type == "LCG"
        oracle = FrankWolfe.lazified_conditional_gradient
    elseif oracle_type == "BCG"
        oracle = FrankWolfe.blended_conditional_gradient
    elseif oracle_type == "BPCG"
        oracle = FrankWolfe.blended_pairwise_conditional_gradient
    end
    
    # create L1 ball as feasible region
    region = FrankWolfe.LpNormLMO{1}(tau-1)


    # call oracles
    if inverse_hessian_boost in ["weak", "full"]
        display("Inverse Hessian Boosting (IHB) is active. Vanilla Frank-Wolfe is used for the IHB run.")

        # compute starting point for IHB
        x0 = l1_projection(solution; radius=tau-1)
        x0 = reshape(x0, length(x0))

        # IHB oracle call
        coefficient_vector, _ = FrankWolfe.frank_wolfe(f, grad!, region, x0; epsilon=epsilon, max_iteration=max_iters)
        if typeof(coefficient_vector) <: FrankWolfe.ScaledHotVector
            coefficient_vector = convert(Vector, coefficient_vector)
        end
        coefficient_vector = vcat(coefficient_vector, [1])
    
        loss = 1/m * norm(data_with_labels * coefficient_vector, 2)^2

        # attempt to find sparse solution if IHB solution found
        if inverse_hessian_boost == "weak" && loss <= psi 
            display("IHB solution found. Attempting to find sparse solution.")

            x0 = compute_extreme_point(region, zeros(Float64, n))
            x0 = Vector(x0)

            tmp_coefficient_vector, _ = oracle(f, grad!, region, x0; epsilon=epsilon, max_iteration=max_iters)
            tmp_coefficient_vector = vcat(tmp_coefficient_vector, [1])
            
            loss2 = 1/m * norm(data_with_labels * tmp_coefficient_vector, 2)^2
            
            if loss2 <= psi
                loss = loss2
                coefficient_vector = tmp_coefficient_vector
            end
        end
    else
        # compute starting vertex
        x0 = compute_extreme_point(region, zeros(Float64, n))
        x0 = Vector(x0)

        # oracle call
        coefficient_vector, _ = oracle(f, grad!, region, x0; epsilon=epsilon, max_iteration=max_iters)
        coefficient_vector = vcat(coefficient_vector, [1])
       
        loss = 1/m * norm(data_with_labels * coefficient_vector, 2)^2
    end      
    return coefficient_vector, loss
end


"""
Runs ABM algorithm to find coefficient vector and computes loss.

# Arguments
- 'oracle_type::String': string denoting which oracle to construct
- 'data::Matrix{Float64}': data (O_evaluations)
- 'labels::Vector{Float64}': labels (term_evaluated)
- 'lambda::Union{Float64, Int64}': regularization parameter (if applicable)
- 'data_squared::Matrix{Float64}': squared data 
- 'data_labels::Vector{Float64}': data' * labels 
- 'labels_squared::Float64': labels' * labels
- 'data_squared_inverse::Union{Matrix{Float64}, Nothing}': inverse of data_squared (default is nothing)

# Returns
- 'coefficient_vector::Vector{Float64}': coefficient vector minimizing ABM optimization problem
- 'loss::Float64': loss w.r.t. 'coefficient_vector' 
"""
function abm(
        data::Matrix{Float64}, 
        labels::Vector{Float64},
        data_squared::Matrix{Float64}, 
        data_labels::Vector{Float64},
        labels_squared::Float64
        )    
    data_with_labels = hcat(data, labels)
    m, n = size(data_with_labels)
    
    # prepare data for SVD and perform SVD
    if m > n
        data_squared_with_labels = hcat(data_squared, data_labels)
        bottom_row = vcat(data_labels, labels_squared)
        bottom_row = transpose(bottom_row)
        data_squared_with_labels = vcat(data_squared_with_labels, bottom_row)
        F = svd(data_squared_with_labels)
    else
        F = svd(data_with_labels; full=true)
    end
    # extract coefficient vector
    Vt = F.Vt
    coefficient_vector = Vt[:, end]
    loss = 1/size(data, 1) * norm(data_with_labels * coefficient_vector, 2)^2
    
    return coefficient_vector, loss
end
