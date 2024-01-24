```@meta
EditURL = "../../examples/docs_custom_oracle.jl"
```

# Custom Oracles
In this section we will present an example of how to construct your own oracle to run in $\texttt{OAVI}$. We will be showing this for the objective function $\frac{1}{m}\|Ax + b\|_2^2$ and a Frank-Wolfe oracle. It is easily extendable to other setups.

## Preparing data
Since custom oracles may vary in the data needed for the computations, we only require a custom oracle to need `data` and `labels` to construct everything, i.e. the data matrix $A$ and the label vector $b$. The first step is to prepare the data you need for your objective function. For our example, we extend our function to better see which combinations of the inputs we need.
```math
\frac{1}{m}\|Ax+b\|_2^2 = \frac{1}{m} \langle Ax+b, Ax+b\rangle = \frac{1}{m} (\langle Ax, Ax\rangle + 2\langle Ax, b\rangle + \langle b, b\rangle).
```
Writing this in the form of matrix-vector multiplication, we get
```math
\frac{1}{m} (\langle Ax, Ax\rangle + 2\langle Ax, b\rangle + \langle b, b\rangle) = \frac{1}{m}(x^\top A^\top A x + 2x^\top A^\top b + b^\top b).
```
Let's write a function that computes all the necessary parts. Inlcuding a $\frac{2}{m}$ factor instead of $\frac{1}{m}$ is just personal preference. We simply divide by $2$ later.

````@example docs_custom_oracle
using AVI
using FrankWolfe
using LinearAlgebra
using Random

function prepare_data(A, b)
    m, _ = size(A)
    # A.T A
    A_squared = 2/m * (A' * A)

    # A.T b
    A_b = 2/m * (A' * b)

    # b.T b
    b_squared = 2/m * (b' * b)

    return A_squared, A_b, b_squared
end
````

## Creating the objective function
The next step is to create the objective function (and in our case its gradient). For this we use the previously defined function to prepare the data and create the objective.

````@example docs_custom_oracle
function objective(data, labels)
    # prepare data for objective
    A_squared, A_b, b_squared = prepare_data(data, labels)

    # define objective
    function evaluate_objective(x)
        return (1 / 2) * (x' * A_squared * x) + (x' * A_b) + (1 / 2) * b_squared
    end

    # define gradient
    function evaluate_gradient!(storage, x)
        # 'storage' is specific to structure in 'FrankWolfe.jl'
        return storage .= A_squared * x + A_b
    end

    return evaluate_objective, evaluate_gradient!
end
````

## Defining the feasible region
Since we are dealing with a constrained optimization problem, defining a feasible region is necessary. We will use FrankWolfe.jl's built-in feasible regions, however one may construct their own feasible region through other means, so for the sake of this presentation we will once again define a function. This function will define the $\ell_p$-ball with a given radius.

````@example docs_custom_oracle
function feasible_region(p, radius)
    return FrankWolfe.LpNormLMO{p}(radius)
end
````

Closely tied with the feasible region is the choice of a starting point for the oracle. For the `LpNormLMO` we can use `compute_extreme_point` to obtain a starting vertex of our feasible region by calling `region = feasible_region(1, 1000)` and after that `x0 = compute_extreme_point(region, zeros(Float64, size(data, 2)))`

## Calling the Oracle
The final step for our custom oracle is calling the oracle with all the things we prepared. $\texttt{OAVI}$ requires degree-lexicographical term ordering and assumes a leading term coefficient of $1$. Hence, your oracle should find a coefficient vector $x$ that only contains coefficients for `data` as the coefficient for `labels` is fixed at $1$. With the coefficient vector obtained, we compute the loss and return both `coefficient_vector` and `loss`. As an example we take the `blended_conditional_gradient` algorithm as an oracle.

````@example docs_custom_oracle
function custom_oracle(data, labels; epsilon=1.0e-7, max_iteration=10000)
    # get necessary data
    f, grad! = objective(data, labels)
    region = feasible_region(1, 1000.)
    x0 = compute_extreme_point(region, zeros(Float64, size(data, 2)))

    # obtain coefficient vector
    coefficient_vector, _ = blended_conditional_gradient(
                                                        f,
                                                        grad!,
                                                        region,
                                                        x0;
                                                        epsilon=epsilon,
                                                        max_iteration=max_iteration
                                                        )

    # compute loss
    loss = f(coefficient_vector)

    # extend coefficient vector with leading term coefficient 1
    coefficient_vector = vcat(coefficient_vector, [1])

    return coefficient_vector, loss
end
````

## Putting it all together
Now that you know what your custom oracle is expected to return and how one can go about defining such a constructor, it remains to give it to $\texttt{OAVI}$ to use.

````@example docs_custom_oracle
# some data
X = rand(10000, 10)

# transformation obtained through OAVI by using your custom oracle
X_transformed, sets = fit_oavi(X; oracle=custom_oracle)
````

You can also pass further kwargs used by your oracle and we will pass them through to your oracle call. This is done by the `oracle_kwargs` argument. Let's say we want to pass the keyword arguments `epsilon=1.0e-5` and `max_iteration=5000` along to our custom oracle.

````@example docs_custom_oracle
kwargs = [(:epsilon, 1.0e-5), (:max_iteration, 5000)]

X_transformed, sets = fit_oavi(X; oracle=custom_oracle, oracle_kwargs=kwargs)
````

The above code will call `custom_oracle` with the keyword arguments `epsilon` and `max_iteration` exchanged for `1.0e-5` and `5000`, respectively, that is, `custom_oracle(data, labels; epsilon=1.0e-5, max_iteration=5000)` and ultimately when calling the algorithm in our example `blended_conditional_gradient(f, grad!, region, x0; epsilon=1.0e-5, max_iteration=5000)`.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

