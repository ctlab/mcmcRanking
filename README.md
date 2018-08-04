
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcmcRanking

[![Travis-CI Build
Status](https://travis-ci.org/ctlab/mcmcRanking.svg?branch=master)](https://travis-ci.org/ctlab/mcmcRanking)
[![Coverage
status](https://codecov.io/gh/ctlab/mcmcRanking/branch/master/graph/badge.svg)](https://codecov.io/github/ctlab/mcmcRanking?branch=master)

## Overview

Tool for estimate probabilities of vertices being in active module using
its likelihoods and it proposes methods for ranking vertices in order of
importance. Estimating probabilities based on Markov chain Monte Carlo
(MCMC) methods.

## Installation

You can install mcmcRanking from github with:

``` r
# install.packages("devtools")
devtools::install_github("ctlab/mcmcRanking")
```

## Illustration

An animation of the MCMC algorithm performing. Vertices are colored
depending on its likelihood, burgundy and gray colors corresponds to
high and low likelihoods respectively. Yellow subgraph is an active
module.

<a href="https://gist.github.com/javlon/36d0ade4be626a3fbb8cc34ee1ac1d69"><img src="https://gist.githubusercontent.com/javlon/36d0ade4be626a3fbb8cc34ee1ac1d69/raw/836b53ae53e1d7f61cdc2bd119f0698c8a50ddf1/mcmc_sample.gif" width="630" height="630"/></a>
