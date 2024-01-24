# How does it work?

`AVI.jl` contains routines to construct approximately vanishing polynomials from the data processed previously. Through the iterations, the data used for construction is expanded with the evaluation of non-leading terms over the data points $\texttt{X}$. 

## Oracle Approximate Vanishing Ideal algorithm 

The Oracle Approximate Vanishing Ideal ($\texttt{OAVI}$) algorithm solves the convex optimization problem 
```math
x^* \in \argmin_{\|x\|_1 \le \tau} \|Ax + b\|_2^2,
```
where $A$ is the matrix of non-leading terms evaluated over the data $\texttt{X}$ and $b$ is the current leading term candidate. Through the use of an oracle, the algorithm finds the coefficient vecotr $x^{*}$ which minimizes the above problem. 

The algorithm is implemented in the [`fit_oavi`](@ref) function. See [E. Wirth and S. Pokutta (2022)](https://proceedings.mlr.press/v151/wirth22a.html) for more details about the method.

### Frank-Wolfe algorithms as Oracle

One such class of methods are 'Frank-Wolfe' algorithms. The [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl/tree/master) package provides implementations of many different version of the Frank-Wolfe algorithm, as well as implementations of Linear Minimization Oracles that find an optimal vertex of the feasible set along a given direction.

You can choose from the standard `frank_wolfe`, the `blended_conditional_gradient` and the `blended_pairwise_conditional_gradient` Frank-Wolfe algorithms by using the keyword argument `oracle = "CG"`, `oracle = "BCG"` and `oracle = "BPCG"`, respectively. For all of those we create the `L1-ball` of radius $\tau$ as well as the squared euclidean norm $\|Ax+b\|_2^2$ as the objective function, as stated above.

### Approximate Buchberger-Möller ($\texttt{ABM}$)

$\texttt{ABM}$ is another related algorithm which shares the degree-lexicographical term ordering as well as processing terms of a given degree $d$ one by one instead of all at once. Due to these similarities we elected to implement $\texttt{ABM}$ inside of [`fit_oavi`](@ref). The algorithm can be run by providing the keyword argument `oracle = "ABM"`. See [J. Limbeck (2013), Chapter 4](https://www.researchgate.net/publication/283651363_Computation_of_Approximate_Border_Bases_and_Applications) for more details about the method.

### Vanishing Component Analysis ($\texttt{VCA}$)
$\texttt{VCA}$ is the final algorithm implemented in this repository. It is once again a method to obtain a polynomial feature transformation with small but significant differences compared to $\texttt{OAVI}$. The main differences are $\texttt{VCA}$ not relying on degree-lexicographic ordering of terms when processing terms and using SVD to obtain coefficient vectors in contrast to the use of an oracle in $\texttt{OAVI}$. See [Livni et al. (2013)](https://proceedings.mlr.press/v28/livni13.html) for additional information.
