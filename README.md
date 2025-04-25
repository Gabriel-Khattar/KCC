
# Quantifying Community Keystoneness with the KCC package

<!-- badges: start -->
<!-- badges: end -->
We introduce a new framework that leverages observational data on community composition to quantify community keystoneness, i.e., the role of a local community in maintaining metacommunity structure over time. 

This framework also offers a valuable tool for practitioners to assess the value of local communities and their contributions to metacommunity dynamics. 

This package contains all functions necessary to estimate community keystoneness using compositional data (i.e., species X sites X time matrices), including its main function `Keystone_calculation`


For practical examples of how to apply this framework, refer to Section 4 of Khattar and Peres‑Neto (2025) “Quantifying Community Keystoneness in Metacommunities under Disturbance” (Oikos, https://doi.org/10.1002/oik.11163) and its Supplementary Information.

## Installation

You can install the development version of KCC with:

``` r
# install devtools if you don’t have it
install.packages("devtools")

# then install the package
devtools::install_github("Gabriel-Khattar/KCC")

# Operational Example
library(vegan)
data (mite)

Keystoneness_calculation(x=mite,rep=10,ranking = T,ncores = 1)


```




