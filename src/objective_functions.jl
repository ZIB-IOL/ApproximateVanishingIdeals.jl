"""
Creates and returns objective function, gradient and solution (through inversion) w.r.t. L2 loss.

# Arguments
- 'data::Matrix{Float64}': evaluations of O_terms 
- 'labels::Vector{Float64}': current border term evaluated
- 'lambda::Union{Float64, Int64}': regularization parameter
- 'data_squared::Matrix{Float64}': data' * data 
- 'data_labels::Vector{Float64}': data' * labels
- 'labels_squared::Float64': labels' * labels
- 'data_squared_inverse::Union{Matrix{Float64}, Nothing}': inverse of data_squared for IHB, optional (default is nothing)

# Returns
- 'solution::Vector{Float64}': solution to (unconstrained) minimization problem
- 'evaluate_function<:Function': objective function
- 'evaluate_gradient!<:Function': gradient of objective function, adjusted to FrankWolfe requirements
"""
function L2Loss(
        data::Matrix{Float64}, 
        labels::Vector{Float64},
        lambda::Union{Float64, Int64}, 
        data_squared::Matrix{Float64}, 
        data_labels::Vector{Float64},
        labels_squared::Float64;
        data_squared_inverse::Union{Matrix{Float64}, Nothing}=nothing
    )

    A = data
    m, n = size(data)
    b = labels
        
    # A_squared
    A_squared = 2/m * data_squared

    if lambda != 0.
        A_squared = A_squared + lambda * Matrix(I, n, n)
    end
        
    # A_b
    A_b = 2/m * data_labels

    # A_squared_inv
    A_squared_inv = nothing
    solution = nothing
    if data_squared_inverse !== nothing
        A_squared_inv = m/2 * data_squared_inverse
        solution = - data_squared_inverse * data_labels
        @assert lambda == 0. "Regularization not implemented for hessian-based algorithms."
    end
        
    # b_squared
    b_squared = 2/m * labels_squared
        
    function evaluate_function(x)
        return ((1 / 2) * (x' * A_squared * x) .+ (A_b' * x) .+ (1 / 2) * b_squared)
    end

    function evaluate_gradient!(storage, x)
        return storage .= A_squared * x + A_b
    end
    
    return solution, evaluate_function, evaluate_gradient!
end
