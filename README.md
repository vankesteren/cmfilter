<p align="center">
  <img src="cmfilter.png" width="300px"></img>
  <br/>
  <span>
    <a href="https://CRAN.R-project.org/package=cmfilter"><img src="http://www.r-pkg.org/badges/version/cmfilter"></img></a>
    <a href="https://travis-ci.org/vankesteren/cmfilter"><img src="https://travis-ci.org/vankesteren/cmfilter.svg?branch=master"></img></a>
    <a href="https://ci.appveyor.com/project/vankesteren/cmfilter"><img src="https://ci.appveyor.com/api/projects/status/f0hbgmqlgkqhdstj?svg=true"></img></a>
  </span>
  <h5 align="center">Coordinate-wise Mediation Filter</h5>
</p>
<br/>

## Description
An `R` package for simultaneous discovery of multiple mediators in an _x → M → y_ system using Coordinate-wise Mediation Filtering.

Keywords: `high-dimensional data`, `feature selection`, `structural equation modeling`, `mediation analysis`

## Installation
The package can be installed from the `r-universe`:

```r
# Enable repository from vankesteren
options(repos = c(
  vankesteren = "https://vankesteren.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Download and install cmfilter in R
install.packages("cmfilter")
```


To install the development version of `cmfilter`, run

```r
devtools::install_github("vankesteren/cmfilter@devel")
```

## Usage
The built-in documentation (run `help(cmf)`) gives information on how to use this package.

## Citation

```
van Kesteren, E. J., & Oberski, D. L. (2019). Exploratory mediation analysis with many potential mediators. Structural Equation Modeling: A Multidisciplinary Journal, 26(5), 710-723.
```
