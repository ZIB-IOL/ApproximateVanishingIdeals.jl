# Obtaining a transformation
With basics about the method out of the way, let's take a look at how we can obtain a feature transformation using the provided algorithms.

## $\texttt{OAVI}$ transformation
To obtain a simple transformation, one can proceed to follow the steps presented earlier. As explained, however, we recommend tweaking a couple of parameters to achieve better results. The two main arguments this applies to are $\psi$, the vanishing parameter, under the keyword `psi` and $\varepsilon$, the accuracy to which the oracle solves the problem, under the keyword `epsilon`. These are by default `psi=0.1` and `epsilon=1.0e-7`. Let's say we want to use `psi=0.005` and `epsilon=1.0e-5`, then the code should look as follows:
```julia
using AVI

X = your_data_set
transformed_data, sets = fit_oavi(X; psi=0.005, epsilon=1.0e-5)
```
While `epsilon` should be adjusted according to the solver used, or personal preference for that matter, the `psi` keyword has a direct influence on the degree of polynomials constructed. It was shown in [Wirth, Kera and Pokutta](https://openreview.net/forum?id=3ZPESALKXO) that the algorithm terminates after having constructed polynomials of degree
```math
D = \lceil -\log(\psi)/\log(4) \rceil,
```
if $\tau \ge (3/2)^D$. The `tau` parameter by default takes the value of this lower bound, but it is possible to adjust it if one wants to. 

### Different built-in Oracles
We previously alluded to three different oracles being available in the code, `frank_wolfe`, `blended_conditional_gradient` and `blended_pairwise_conditional_gradient`. The default option is `frank_wolfe`, which has the keyword arugment `oracle="CG"`. The others can be accessed as follows:
```julia
# Standard Frank-Wolfe. Both these formulations are the same.
fit_oavi(X)
fit_oavi(X; oracle="CG")   

# Blended conditional gradients
fit_oavi(X; oracle="BCG")

# Blended pairwise conditional gradients
fit_oavi(X; oracle="BPCG")
```
### Custom Oracles
It is also possible to provide your own oracle for solving the convex optimization problem. The code expects your constructor to work with the two variable arguments `data` and `labels`, which are the data matrix $A$ and the label vector $b$ found in $\|Ax + b\|^2$. Further, you can add your keyword arguments under `oracle_kwargs`. In the `Examples` section, we will show you how to construct a custom oracle by showing an example with an oracle from [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl/tree/master).

### $\texttt{ABM}$ 
Through the keyword `oracle` we can also call the $\texttt{ABM}$ algorithm instead. This can be accessed as follows:
```julia
# Approximate Buchberger-Möller
fit_oavi(X; oracle="ABM")
```

## $\texttt{VCA}$ transformation
Apart from $\texttt{OAVI}$ with its multiple different possibilities of running, we also provide an implementation of the $\texttt{VCA}$ algorithm introduced in [Livni et al. (2013)](https://proceedings.mlr.press/v28/livni13.html). The algorithm is implemented in [`fit_vca`](@ref) and a feature transformation based on $\texttt{VCA}$ can be obtained by calling
```julia
using AVI

X = your_data_set
X_transformed_vca, sets_vca = fit_vca(X)
```
Similar to $\texttt{OAVI}$, the vanishing parameter $\psi$ is adjustable for $\texttt{VCA}$ through the keyword `psi`. To the best of our knowledge, there is currently no similar guarantee on the degree of constructed polynomials for $\texttt{VCA}$ as there is for $\texttt{OAVI}$, but adjusting `psi` and comparing results is nonetheless recommended.
