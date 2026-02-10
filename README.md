# mixedGP

Code for "Errors-in-variables Gaussian processes for mixed-input regression".

The "main-code" folder contains the R scripts for main functions.

* basics-nominal.R: functions for one nominal input.
* basics-ordinal.R: functions for one ordinal input.
* basics-common.R: functions that are used for both nominal and ordinal inputs.
* plot-config.R: configuration for the plots in the R Markdown files.



The "markdown" folder contains the R Markdown files for conducting the numerical experiments.

* simu-ordinal.Rmd \& simu-nominal.Rmd: synthetic data generation, hyperparameter estimation, and implementation of the Gibbs sampler for latent variables for the ordinal and nominal inputs, respectively.
* mixing-ordinal.Rmd \& mixing-nominal.Rmd: mixing diagnosis, model fitting and prediction performance for the ordinal and nominal inputs, respectively.
