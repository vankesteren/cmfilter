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

<h6 align="center"> Package <code>cmfilter</code> is work in progress. Check back later for more info. </h6>
<br/>



## Description
An `R` package for simultaneous discovery of multiple mediators in an _x → M → y_ system using Coordinate-wise Mediation Filtering.

Keywords: `high-dimensional data`, `feature selection`, `structural equation modeling`, `mediation analysis`

## Installation
The package is not yet available on `CRAN`. To install the package directly from this repository, install the `devtools` package, make sure you have [`R Build Tools` (Windows)](https://cran.r-project.org/bin/windows/Rtools/) installed, and then run the following command:
```r
devtools::install_github("vankesteren/cmfilter")
```

If you have only installed the toolchain for your current architecture (32-bit or 64-bit only), run `options(devtools.install.args = "--no-multiarch")` before installing.


To install the development version of `cmfilter`, run

```r
devtools::install_github("vankesteren/cmfilter@devel")
```

## Usage
The built-in documentation (run `help(cmf)`) gives information on how to use this package. More extensive documentation is under development.
