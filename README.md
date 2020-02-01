## samTEMsel (R package version 0.1.0)
Sparse Additive models for Treatment Effect-Modifier Selection 

An implementation of *a constrained sparse additive regression* for modeling interaction effects between a categorical treatment variable and a set of pretreatment covariates on a scalar-valued outcome; the regression simultaneousely conducts treatment effect-modifier variable selection. The method can effectively identify treatment effect-modifiers exhibiting possibly nonlinear interactions with the treatment. The selected pretreatment characteristics and the associated nonzero component functions can be used as a new set of data-driven features for making individualized treatment recommendations in further analysis. We refer to Park, Petkova, Tarpey, and Ogden (2020) <doi:10.1016/j.jspi.2019.05.008> and Park, Petkova, Tarpey, and Ogden (2020) "A constrained sparse additive model for treatment effect-modifier selection" (pre-print) for detail of the method. The wrapper function of this package is cv.samTEMsel().


#### Description

* samTEMsel - `samTEMsel` main function
* cv.samTEMsel - `samTEMsel` cross-validation function for tuning parameter selection 
* predict_samTEMsel - `samTEMsel` prediction function
* make_ITR - make individualized treatment recommendations (ITRs) based on a `samTEMsel` object
* plot_samTEMsel -  plot component functions from a `samTEMsel` object 


#### To run: 
To install an R package, start by installing the "devtools" package (from CRAN). On R, type: 
```
install.packages("devtools")  # install the devtools package from CRAN
library(devtools)
```

To install the "samTEMsel" package from github, type: 
```
devtools::install_github("syhyunpark/samTEMsel")  # install the samTEMsel package from github 
library(samTEMsel)   # load the samTEMsel package to R 
```

To see some of the example codes appearing in the "help" menu, type:  
```
?samTEMsel   
?cv.samTEMsel
```

