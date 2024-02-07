# How does it work?

`AVI.jl` contains routines to construct approximately vanishing polynomials from data processed previously. Through the iterations, the data used for construction is expanded with the evaluation of non-leading terms over the data points $\texttt{X}$. 

## Oracle Approximate Vanishing Ideal algorithm ($\texttt{OAVI}$)

The Oracle Approximate Vanishing Ideal ($\texttt{OAVI}$) algorithm solves the convex optimization problem 
```math
x^* \in \underset{\|x\|_1 \le \tau}{\textnormal{arg}\min} \|Ax + b\|_2^2,
```
where $A$ is the matrix of non-leading terms evaluated over the data $\texttt{X}$ and $b$ is the current leading term candidate. Through the use of an oracle, the algorithm finds the coefficient vector $x^{*}$ which minimizes the above problem. This $x^{*}$ is then used to determine whether the constructed polynomial is a vanishing one, or whether the leading term candidate does not produce a vanishing polynomial. 

The formulation above may look similar to a regression problem, but that is not what is happening here. The convex optimization problem is solved independently for each distinct leading term candidate and only to minimize the distance to the zero-polynomial for this specific combination of non-leading terms $A$ and leading term $b$. So, for each leading term candidate we wish to find 
```math
x^* \in \underset{\|x\|_1 \le \tau}{\textnormal{arg}\min} \|(Ax + b) - \mathbf{0}\|_2^2,
```
where $Ax+b$ is describes the polynomial with non-leading terms $A$ and leading term $b$ and $\mathbf{0}$ is the zero-polynomial, instead of minimizing, e.g. the average distance of our polynomials to the zero-polynomial 

The algorithm is implemented in the [`fit_oavi`](@ref) function. See [E. Wirth and S. Pokutta (2022)](https://proceedings.mlr.press/v151/wirth22a.html) for more details about the method.

### Frank-Wolfe algorithms as Oracle

One such class of oracles that fit the optimization problem in $\texttt{OAVI}$ are 'Frank-Wolfe' algorithms. The [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl/tree/master) package provides implementations of many different version of the Frank-Wolfe algorithm, as well as implementations of Linear Minimization Oracles that find an optimal vertex of the feasible set along a given direction.

You can choose from the standard `frank_wolfe`, the `blended_conditional_gradient` and the `blended_pairwise_conditional_gradient` Frank-Wolfe algorithms by using the keyword argument `oracle = "CG"`, `oracle = "BCG"` and `oracle = "BPCG"`, respectively. For all of those we create the `L1-ball` of radius $\tau$ as well as the squared euclidean norm $\|Ax+b\|_2^2$ as the objective function, as stated above.

### Approximate Vanishing Ideal ($\texttt{AVI}$)
The basis of $\texttt{OAVI}$ and many related methods is the Approximate Vanishing Ideal ($\texttt{AVI}$) algorithm introduced in [Heldt, Kreuzer and Pokutta (2009)](https://www.sciencedirect.com/science/article/pii/S0747717109000935). $\texttt{AVI}$ itself is based on the (exact) Buchberger-Möller ($\texttt{BM}$) algorithm [Möller and Buchberger (1982)](https://link.springer.com/chapter/10.1007/3-540-11607-9_3) which proved unreliable if the data set was corrupted with noise. To handle this instability, Heldt et al. abstained from computing an exact basis of vanishing polynomials and, by allowing the polynomials to only vanish _approximately_, created a numerically stable algorithm to compute the _approximate_ vanishing ideal. We do not provide an explicit implementation of $\texttt{AVI}$, as more and more modern methods building upon it surfaced, which proved to produce even better numerical stability and at the same time retaining many of the original ideas of $\texttt{AVI}$. This acknowledgement is not intended to diminish the algorithm's significance but instead to highlight the impressive strides made by many different researchers in this field to refine and innovate upon the foundation laid by the $\texttt{AVI}$ algorithm.

### Approximate Buchberger-Möller ($\texttt{ABM}$)

The Approximate Buchberger-Möller ($\texttt{ABM}$) algorithm is, next to $\texttt{OAVI}$, another related algorithm which shares the degree-lexicographical term ordering as well as processing terms of a given degree $d$ one by one instead of all at once. Due to these similarities we elected to implement $\texttt{ABM}$ inside of [`fit_oavi`](@ref). The algorithm can be run by providing the keyword argument `oracle = "ABM"`. See [J. Limbeck (2013), Chapter 4](https://www.researchgate.net/publication/283651363_Computation_of_Approximate_Border_Bases_and_Applications) for more details about the method.

## Vanishing Component Analysis ($\texttt{VCA}$)
The Vanishing Component Analysis ($\texttt{VCA}$) algorithm is the final algorithm implemented in this repository. It is once again a method to obtain a polynomial feature transformation with small but significant differences compared to $\texttt{OAVI}$. The main differences are $\texttt{VCA}$ not relying on degree-lexicographic ordering of terms when processing terms and using SVD to obtain coefficient vectors in contrast to the use of an oracle in $\texttt{OAVI}$. The lack of a degree-lexicographic ordering often leads to the coefficient vectors exploding in magnitude. See [Livni et al. (2013)](https://proceedings.mlr.press/v28/livni13.html) for additional information.
