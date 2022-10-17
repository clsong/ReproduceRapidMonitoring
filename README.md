
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ReproduceRapidMonitoring

<!-- badges: start -->
<!-- badges: end -->

The goal of this package is to reproduce the paper “Rapid monitoring for
ecological persistence” by Chuliang Song, Benno Simmons, Marie-Josée
Fortin, Andrew Gonzalez, Christopher Kaiser-Bunbury, Serguei Saavedra.

## Installation

You can install the development version of ReproduceRapidMonitoring from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clsong/ReproduceRapidMonitoring")
```

## Reproduce

The main figures can be reproduced simply by running the corresponding
functions. For example, figure 3B can be reproduced as

``` r
library(ReproduceRapidMonitoring)
generate_figure3B_plot()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
