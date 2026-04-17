# sage_periods

This repository hosts the implementation of a SageMath package containing
algorithms for computing annihilating linear differential equations
for diagonals and period integrals of multivariate rational functions.

The package works with any reasonably recent version of SageMath, we
recommend to have SageMath 10.1 or newer.
Documentation is available at <https://acsvmath.github.io/sage_periods/>.

Our methods are based on the [MAGMA period package](https://github.com/lairez/periods) 
of [Pierre Lairez](https://mathexp.eu/lairez/) which implements the algorithm described 
in the paper *[Computing Periods of Rational 
Integrals](https://www.ams.org/journals/mcom/2016-85-300/S0025-5718-2015-03054-3/)*, and integrate 
with the [sage-acsv](https://github.com/ACSVMath/sage_acsv) and 
[ore-algebra](https://github.com/mkauers/ore_algebra/) SageMath packages.


## Quickstart

To install the package from the source code, clone the git repository and run the command
```sh
sage -pip install .
```
from the root directory (the directory containing the `pyproject.toml` file).

For development, use `sage -pip install -e .` for an editable installation.

Alternatively, to install the latest version of the main branch directly from
the GitHub repository, run
```sh
sage -pip install git+https://github.com/ACSVMath/sage_periods.git
```

## Example of Use

After installation the computation
```

sage: from sage_periods import compute_diagonal_annihilator
sage: var('x y')
sage: F = 1/(1-x-y)
sage: compute_diagonal_annihilator(F)
(t - 1/4)*Dt + 1/2

``` 
shows that the main diagonal 

$$ S(t) = \sum_{n \geq 0}\binom{2n}{n}t^n = (1-4t)^{-1/2}$$ 

of $1/(1-x-y)$ satisfies $(t-1/4)S'(t) + (1/2)S(t) = 0$.

Further examples can be found in the [documentation](https://acsvmath.github.io/sage_periods/) for each function.

