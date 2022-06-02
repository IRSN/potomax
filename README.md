
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
distribution, the model include the three parameters
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
(rate of the Poisson process for the exceedances over the threshold
![u](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u "u")),
and the two GP parameters
![\\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma "\sigma")
(scale) and
![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
(shape). The vector of parameters
![\\boldsymbol{\\theta} = \[\\lambda, \\, \\sigma,\\, \\xi\]^\\top](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Ctheta%7D%20%3D%20%5B%5Clambda%2C%20%5C%2C%20%5Csigma%2C%5C%2C%20%5Cxi%5D%5E%5Ctop "\boldsymbol{\theta} = [\lambda, \, \sigma,\, \xi]^\top")
can be transformed into the vector
![\\boldsymbol{\\theta}^\\star = \[\\mu^\\star, \\, \\sigma^\\star, \\,\\xi^\\star\]^\\top](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Ctheta%7D%5E%5Cstar%20%3D%20%5B%5Cmu%5E%5Cstar%2C%20%5C%2C%20%5Csigma%5E%5Cstar%2C%20%5C%2C%5Cxi%5E%5Cstar%5D%5E%5Ctop "\boldsymbol{\theta}^\star = [\mu^\star, \, \sigma^\star, \,\xi^\star]^\top")
of so-called “Poisson-Process” (PP) parameters that describe the GEV
distribution of the maximum on a block with a fixed reference duration
![w^\\star](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%5E%5Cstar "w^\star"),
usually taken as one year. The PP parameters do not depend on the
threshold
![u](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u "u")
but they depend on
![w^\\star](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%5E%5Cstar "w^\star").
The Poisson-GP parameters do not depend on
![w^\\star](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%5E%5Cstar "w^\star")
but depend on the threshold
![u](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u "u").
The shape parameters
![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
and
![\\xi^\\star](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi%5E%5Cstar "\xi^\star")
are identical in the two parameterisations. It is often simpler to use
the Poisson-GP parameterisation because the rate
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
can be concentrated out of the likelihood.

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

## Version 0.2.1

-   New vignette: *R Package potomax: Overview*, including the Venice
    example.

## Version 0.2.0

-   The one-parameter exponential distribution can now be used as a
    distribution for the excesses. This leads to Gumbel-distributed
    block maxima.

-   Bug fix. The quantile function `qGPD2` did not work as expected with
    `lower.tail = FALSE`.

-   The approximation of the GPD2 distribution functions for
    ![\\xi \\approx 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi%20%5Capprox%200 "\xi \approx 0")
    has been improved. For each of the three functions `dGPD2`, `pGPD2`
    and `qGPD2`, a quadratic function of
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
    is used for small
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
    providing a Taylor series approximation. The derivatives w.r.t. the
    parameters are consistently taken as functions of
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
    that are linear (gradient) and constant (hessian). The distribution
    functions are now smooth functions of
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")
    even for very small
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi"),
    say
    ![\\xi \\approx 10^{-6}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi%20%5Capprox%2010%5E%7B-6%7D "\xi \approx 10^{-6}")
    or less. Some tests devoted to this “small
    ![\\xi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi "\xi")”
    approximation have been added to the package.

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
