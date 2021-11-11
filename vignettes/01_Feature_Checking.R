#' ---
#' title: "Feature Checking (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#+ setup, include=FALSE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varchk <- "out/vig1/"
file_varchk <- file.path( dir_varchk, "tcga_ucec_cna_cnv_varchk.rds" )

# load data
load( "data/tcga_ucec_cna_cnv.Rdata" )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars

# Run variable checks
start_time <- Sys.time()
vck <- check_variables( dat, re, va )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Variable checking ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )

# Save results and render report
saveRDS( vck, file = file_varchk )
render_variable_checks_summary( vck, dir_varchk )


#' # Introduction
#' 
#' Building predictive models requires as an immediate first step the inpsection of available data:
#' \describe{
#' \item{Missing values}{It's key to understand the level of *missing values* per variable. Multi-variate
#' models may fail to train as the number of complete cases may drop. Alternatively, approaches for
#' handling missing data may be considered. If key variables bear too many missing values, this may
#' be indicative of a relevant process or also of data quality.}
#' \item{Data type and Information Content}{Regardless of data type, a variable should comprise a certain level
#' of information, that is should attain more than one value or not be constant except for a few cases. This
#' property can be estimated as the variable *Entropy* and the *Mutual Information* with respect to a response.}
#' \item{Outliers and Extreme Cases}{For numerical variables, it is essential to understand the propability and
#' extend of potential *outliers*. Some modelling approaches may be severely impacted by those. We can use two
#' approaches: 1) is based on the so-called _Robust Outlier Test_ using the range outsider [Q1 - c * IQR, Q3 + c * IQR]} 
#' where c = 1.5 or 3.0 to define outliers, and 2) a statistical outlier test like the Chi-squared outlier test, Dixon or 
#' Grubbs test providing p-values for the likelihood that one tail comprises an outlier. For the Gibbs test this is based
#' on the ratio of sample range divided by sample standard deviation (Gibbs U-Statistic).
#' }
#' 


