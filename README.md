# samTEMsel
Sparse Additive models for Treatment Effect-Modifier Selection 

An implementation of a constrained sparse additive regression for modeling interaction effects between a categorical treatment variable and a set of pretreatment covariates on a scalar-valued outcome; the regression simultaneousely conducts treatment effect-modifier variable selection. The method can effectively identify treatment effect-modifiers exhibiting possibly nonlinear interactions with the treatment. The selected pretreatment characteristics and the associated nonzero component functions can be used as a new set of data-driven features for constructing ITRs in further analysis. We refer to Park, Petkova, Tarpey, and Ogden (2020) <doi:10.1016/j.jspi.2019.05.008> and Park, Petkova, Tarpey, and Ogden (2020) "A constrained sparse additive model for treatment effect-modifier selection" (pre-print) for detail of the method. The wrapper function of this package is cv.samTEMsel().


## Descrpition of R functions 


* cv.samTEMsel - `samTEMsel` cross-validation function
* samTEMsel - `samTEMsel` main function
* predict_samTEMsel - `samTEMsel` prediction function
* make_ITR - make individualized treatment recommendations (ITRs) based on a `samTEMsel` object
* plot_samTEMsel -  plot component functions from a `samTEMsel` object 


## Example R Commands

To install an R package, first need to install the "devtools" package (from CRAN). On the R console, type: 
`
install.packages("devtools")  # install the devtools package from CRAN
library(devtools)
`

To install the "samTEMsel" package to R from Github
`
devtools::install_github("syhyunpark/samTEMsel")  #install the package from Github 
`
Load the samTEMsel package into R, and type ?samTEMsel or ?cv.samTEMsel to see some example codes appearing in the "help" menu 
`
library(samTEMsel)
?samTEMsel
?cv.samTEMsel
`

