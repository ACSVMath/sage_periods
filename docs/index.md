---
icon: lucide/rocket
---

# sage-periods: A SageMath Package to Compute Rational Diagonals and Periods 

This package computes annihilating linear differential equations for diagonals of rational functions and period integrals of rational functions with a single parameter.

Our methods are based on the [MAGMA period package](https://github.com/lairez/periods) of [Pierre Lairez](https://mathexp.eu/lairez/) which implements the algorithm described in the paper *[Computing Periods of Rational Integrals](https://www.ams.org/journals/mcom/2016-85-300/S0025-5718-2015-03054-3/)*, and integrate with the [sage-acsv](https://github.com/ACSVMath/sage_acsv) and [ore-algebra](https://github.com/mkauers/ore_algebra/) SageMath packages.

Most users will need only the functions:

- [`compute_diagonal_annihilator`](reference/#sage_periods.picard_fuchs.compute_diagonal_annihilator) - Takes a symbolic rational function $F(x_1,\dots,x_d)$ and a *direction vector* $\mathbf{r}\in\mathbb{N}^d$ and computes an operator representing a linear ODE with polynomial coefficients annihilating the *$\mathbf{r}$-diagonal* $S(t) = \sum_{n \geq 0}f_{n\mathbf{r}}t^n$ defined by the coefficients $f_{\mathbf{i}}$ of *any* convergent Laurent expansion of $F$. By default the *main diagonal* $\mathbf{r}=\mathbf{1}$ is chosen.

- [`compute_period_annihilator`](reference/#sage_periods.picard_fuchs.compute_period_annihilator) - Takes a symbolic rational function $R(x_1,\dots,x_d,t)$ and computes an operator representing a linear ODE with polynomial coefficients annihilating all *period integrals* $P(t) = \int_\Gamma R(\mathbf{x},t) d\mathbf{x}$ for suitable closed chains of integration $\Gamma$.  

All public facing functions and classes in our package are documented in our [reference manual](reference.md).

## Quickstart

To use the package, simply download the source code from the [GitHub repository](https://github.com/ACSVMath/sage_periods), then import the commands from the package.

```

from sage_periods import compute_diagonal_annihilator

``` 

!!! warning
   
      Make sure you're running SageMath from the same directory where the source code folder is located.

<!-- TODO: Add to Pip. Once we're ready for the 0.1.0 version + have workflows set up for this. -->
<!-- The easiest way to install the latest released version of the package
is via PyPI simply by running

```
sage -pip install sage-periods
```
or, alternatively, executing a cell containing
```
%pip install sage-periods
```
in a SageMath Jupyter notebook. -->

## Examples

After importing `compute_diagonal_annihilator` as above, one can run
```

sage: var('x y')
sage: F = 1/(1-x-y)
sage: compute_diagonal_annihilator(F)
(t - 1/4)*Dt + 1/2

``` 
to show that the main diagonal 

$$ S(t) = \sum_{n \geq 0}\binom{2n}{n}t^n = (1-4t)^{-1/2}$$ 

of $1/(1-x-y)$ satisfies $(t-1/4)S'(t) + (1/2)S(t) = 0$.

Further examples can be found in the documentation for each function.
