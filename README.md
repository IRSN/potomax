
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package potomax

**potomax** is a R package financed by IRSN Behrig.

**potomax** is devoted to Extreme Value analysis and merges the two
classical approaches *POT* (Peaks Over Threshold) and *block MAXima*. In
both cases we may consider that the observations are related in some way
to a *Poisson-GP marked process* i.e., a marked Poisson process with the
marks following the Generalised Pareto (GP) distribution. While this
framework is classical for POT, it is also convenient for block maxima
with arbitrary block durations. This approach allows to censor the block
maxima that are too small to be considered as extreme, as is often
needed for one-year blocks. It also allows to cope with heterogeneous
data as met when using historical information.

For the standard case where the excesses follow the two-parameter GP
distribution, the model include the three parameters *λ* (rate of the
Poisson process for the exceedances over the threshold *u*), and the two
GP parameters *σ* (scale) and *ξ* (shape). The vector of parameters
**θ** = \[*λ*, *σ*, *ξ*\]<sup>⊤</sup> can be transformed into the vector
**θ**<sup>⋆</sup> = \[*μ*<sup>⋆</sup>, *σ*<sup>⋆</sup>, *ξ*<sup>⋆</sup>\]<sup>⊤</sup>
of so-called “Poisson-Process” (PP) parameters that describe the GEV
distribution of the maximum on a block with a fixed reference duration
*w*<sup>⋆</sup>, usually taken as one year. The PP parameters do not
depend on the threshold *u* but they depend on *w*<sup>⋆</sup>. The
Poisson-GP parameters do not depend on *w*<sup>⋆</sup> but depend on the
threshold *u*. The shape parameters *ξ* and *ξ*<sup>⋆</sup> are
identical in the two parameterisations. It is often simpler to use the
Poisson-GP parameterisation because the rate *λ* can be concentrated out
of the likelihood.

While the Poisson-GP formulation was already implemented in the
**Renext** package, **potomax** puts emphasis on the use of
*profile-likelihood* inference rather than on the cheapest widespread
*delta method* that is used in **Renext**. The profile-likelihood
confidence intervals on the parameters and return levels are known to
have better coverage rate than those based on the delta method. However
the profile-likelihood inference requires to repeatedly solve
optimisation problems, with possible numerical issues. The derivatives
of the distribution functions w.r.t the parameters are provided to
facilitate the optimisation tasks. Graphical diagnostics are available
to check that the profile-likelihood results are correct.

The package provides classical S3 methods for fitted models such as
`summary`, `coef`, `logLik`, `confint`, … It also provides the methods
`RL` to compute Return Levels and their confidence intervals, and
`autoplot` method to produce “RL plots”.

# News

## Version 0.2.0

-   The one-parameter exponential distribution can now be used as a
    distribution for the excesses. This leads to Gumbel-distributed
    block maxima.

-   Bug fix. The quantile function `qGPD2` did not work as expected with
    `lower.tail = FALSE`.

-   The approximation of the GPD2 distribution functions for *ξ* ≈ 0 has
    been improved. For each of the three functions `dGPD2`, `pGPD2` and
    `qGPD2`, a quadratic function of *ξ* is used for small *ξ* providing
    a Taylor series approximation. The derivatives w.r.t. the parameters
    are consistently taken as functions of *ξ* that are linear
    (gradient) and constant (hessian). The distribution functions are
    now smooth functions of *ξ* even for very small *ξ*, say
    *ξ* ≈ 10<sup>−6</sup> or less. Some tests devoted to this “small
    *ξ*” approximation have been added to the package.

# Install release version from GitHub

## Using the *devtools* package

Note that if you are using Windows, you need to have the
[Rtools](https://cran.r-project.org/bin/windows/Rtools) installed.
Provided that the **devtools** package is installed you can then in an R
session use

``` r
library(devtools)
install_github("yvesdeville/potomax", dependencies = TRUE, auth_token = myToken)
```

where `myToken` stands for *your* token. This should install the package
and make it ready to use.

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`, see the **devtools** package
documentation.

## Clone, build and install

### Cloning the repository

If you do not have yet a local `potomax` repository, use `git clone` to
clone the `potomax` repository

``` bash
git clone https://github.com/yvesdeville/potomax
```

This will create a `potomax` sub-directory of the current directory,
i.e. the directory from which the git command was issued. Of course this
can work only if you have the authorisation to clone.

### Installation on Unix and MacOs systems

With these systems you can install a package from its source. Move to
the parent directory of your cloned repository and use the following
command from a terminal to create a tarball source file

``` bash
R CMD build potomax
```

This will produce a source tarball `potomax_x.y.z` where `x`, `y` and
`z` stand for the major, minor and patch version numbers. Then you can
install from a command line

``` bash
R CMD INSTALL potomax_x.y.z
```

Note that you must also have all the packages required by **potomax**
installed.

If you are using the **RStudio** IDE, you can alternatively use menus.

### Install and precompile for Windows

In order to install the package from its source, you must have a
suitable Windows platform with
[Rtools](https://cran.r-project.org/bin/windows/Rtools) installed. Then
you can proceed as Unix or MacOS users, with a `build` step from command
line.

If you can not (or do not want to) install the **Rtools** you may get a
trusted binary from a friend or colleague next to you.

### Creating a binary precompiled file for Windows

If you have the **Rtools** installed, you can create a binary. Using a
terminal, move if necessary by using `cd` to the directory containing
the source tarball and R command, and then type

``` bash
R CMD INSTALL --build potomax_x.y.z
```

This will create a `.zip` file that can be used on a Windows platform
which may not be equipped with *Rtools*. For instance, with **RStudio**
you can use the menu `Tools/Install Packages` and select
`Install from:`.

### Precompiled binaries

You can make the resulting binary file available to Windows users who do
not have the **Rtools** installed and are allowed to use the package.
Make sure to conform to the non-disclosure agreement by IRSN and contact
Lise Bardet (`lise.bardet` at `irsn.fr`) in case of doubt.

Also make sure that the **major** version number `x` and the minor
version number `y` are the same as those of the R environment where the
package is to be used.
