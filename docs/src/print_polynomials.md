# Printing polynomials
Apart from just obtaining an $\texttt{OAVI}$ transformation for some data, you might also want to take a look at the polynomials themselves. To do that we provide the function `print_polynomials` which takes as input a `SetsOandG` instance `sets`. The function then puts the polynomials into LaTeX format and prints them, line by line, to the command line.
```julia
# some data
X = rand(100, 3)
_, sets = fit_oavi(X; psi=0.005)

# print polynomials
print_polynomials(sets)
```
As an example, the output may look something like this: 
```julia
x_{3}^{2} - 0.1x_2x_3 + 0.09x_1x_3 + 0.07x_{2}^{2} - 0.05x_1x_2 + 0.04x_{1}^{2} - 1.03x_3 - 0.01x_2 - 0.04x_1 + 0.18
x_{1}^{3} + 0.01x_2x_3 - 0.03x_1x_3 - 0.01x_1x_2 - 1.48x_{1}^{2} + 0.01x_3 + 0.6x_1 - 0.06
x_{1}^{2}x_2 + 0.01x_2x_3 - 0.02x_1x_3 - 0.05x_{2}^{2} - 1.02x_1x_2 - 0.44x_{1}^{2} + 0.01x_3 + 0.22x_2 + 0.46x_1 - 0.09
x_1x_{2}^{2} - 0.03x_2x_3 + 0.08x_1x_3 - 0.52x_{2}^{2} - 1.0x_1x_2 - 0.01x_{1}^{2} - 0.02x_3 + 0.51x_2 + 0.13x_1 - 0.07
x_{2}^{3} + 0.07x_2x_3 + 0.03x_1x_3 - 1.41x_{2}^{2} - 0.04x_1x_2 - 0.06x_{1}^{2} - 0.04x_3 + 0.49x_2 + 0.04x_1 - 0.02
x_{1}^{2}x_3 - 0.96x_1x_3 + 0.02x_{2}^{2} - 0.03x_1x_2 - 0.54x_{1}^{2} + 0.15x_3 - 0.01x_2 + 0.53x_1 - 0.09
x_1x_2x_3 - 0.51x_2x_3 - 0.42x_1x_3 - 0.01x_{2}^{2} - 0.44x_1x_2 - 0.01x_{1}^{2} + 0.2x_3 + 0.24x_2 + 0.2x_1 - 0.1
x_{2}^{2}x_3 - 0.94x_2x_3 - 0.04x_1x_3 - 0.44x_{2}^{2} - 0.03x_1x_2 + 0.02x_{1}^{2} + 0.16x_3 + 0.44x_2 + 0.01x_1 - 0.08
x_1x_{3}^{2} - 0.09x_2x_3 - 0.97x_1x_3 + 0.02x_{2}^{2} + 0.06x_{1}^{2} + 0.02x_3 + 0.01x_2 + 0.11x_1
x_2x_{3}^{2} - 1.04x_2x_3 + 0.02x_{2}^{2} + 0.06x_{1}^{2} - 0.01x_3 + 0.18x_2 - 0.05x_1 + 0.01
```
Then, it is as easy as copy and pasting it into LaTeX with your preferred formatting and you can take a look at them.
\[
x_{3}^{2} - 0.1x_2x_3 + 0.09x_1x_3 + 0.07x_{2}^{2} - 0.05x_1x_2 + 0.04x_{1}^{2} - 1.03x_3 - 0.01x_2 - 0.04x_1 + 0.18\\
x_{1}^{3} + 0.01x_2x_3 - 0.03x_1x_3 - 0.01x_1x_2 - 1.48x_{1}^{2} + 0.01x_3 + 0.6x_1 - 0.06\\
x_{1}^{2}x_2 + 0.01x_2x_3 - 0.02x_1x_3 - 0.05x_{2}^{2} - 1.02x_1x_2 - 0.44x_{1}^{2} + 0.01x_3 + 0.22x_2 + 0.46x_1 - 0.09\\
x_1x_{2}^{2} - 0.03x_2x_3 + 0.08x_1x_3 - 0.52x_{2}^{2} - 1.0x_1x_2 - 0.01x_{1}^{2} - 0.02x_3 + 0.51x_2 + 0.13x_1 - 0.07\\
x_{2}^{3} + 0.07x_2x_3 + 0.03x_1x_3 - 1.41x_{2}^{2} - 0.04x_1x_2 - 0.06x_{1}^{2} - 0.04x_3 + 0.49x_2 + 0.04x_1 - 0.02\\
x_{1}^{2}x_3 - 0.96x_1x_3 + 0.02x_{2}^{2} - 0.03x_1x_2 - 0.54x_{1}^{2} + 0.15x_3 - 0.01x_2 + 0.53x_1 - 0.09\\
x_1x_2x_3 - 0.51x_2x_3 - 0.42x_1x_3 - 0.01x_{2}^{2} - 0.44x_1x_2 - 0.01x_{1}^{2} + 0.2x_3 + 0.24x_2 + 0.2x_1 - 0.1\\
x_{2}^{2}x_3 - 0.94x_2x_3 - 0.04x_1x_3 - 0.44x_{2}^{2} - 0.03x_1x_2 + 0.02x_{1}^{2} + 0.16x_3 + 0.44x_2 + 0.01x_1 - 0.08\\
x_1x_{3}^{2} - 0.09x_2x_3 - 0.97x_1x_3 + 0.02x_{2}^{2} + 0.06x_{1}^{2} + 0.02x_3 + 0.01x_2 + 0.11x_1\\
x_2x_{3}^{2} - 1.04x_2x_3 + 0.02x_{2}^{2} + 0.06x_{1}^{2} - 0.01x_3 + 0.18x_2 - 0.05x_1 + 0.01.\\
\]
The function `print_polynomials` has an additional keyword argument `digits` with which you can decide up to how many digits the coefficients should be rounded, e.g. `print_polynomials(sets; digits=4)` to round to $4$ decimal places. There is a choice you have to make here as some coefficients may be smaller than the decimal place you are rounding to but can nonetheless influence the evaluation of the polynomial in a meaningful way. Be aware of that when printing the polynomials and do not fall under the impression that these are the exact polynomials found during $\texttt{OAVI}$. We completely omit terms that have a rounded coefficient of $0.0$, as well as omitting coefficients that have a rounded value of $1.0$. This is done to declutter the printed polynomials, but, as mentioned above, has the side effect that some terms contained in the polynomial may not be printed due to the coefficients being too small.