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

library( numDeriv )
library( pec )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_ucec_cna_cnv.Rdata")

#' # Introduction
#' 
#' Building predictive models requires as an immediate first step the check of
#' available data:
#' \describe{
#' \item{Missing values}{It's key to understand the level of *missing values* per variable. Multi-variate
#' models may fail to train as the number of complete cases may drop. Alternatively, approaches for
#' handling missing data may be considered. If key variables bear too many missing values, this may
#' be indicative of a relevant process or also of data quality.}
#' \item{Data type and Information Content}{Regardless of data type, a variable should comprise a certain level
#' of information, that is should attain more than one value or not be constant except for a few cases. This
#' property can be estimated as the variable *Entropy*.}
#' \item{Outliers and Extreme Cases}{For numerical variables, it is essential to understand the propability and
#' extend of potential *outliers*. Some modelling approaches may be severely impacted by those. We can use two
#' approaches: 1) is based on the so-called _Robust Outlier Test_ using the range outsider [Q1 - c * IQR, Q3 + c * IQR]} 
#' where c = 1.5 or 3.0 to define outliers, and 2) a statistical outlier test like the Chi-squared outlier test, Dixon or 
#' Grubbs test providing p-values for the likelihood that one tail comprises an outlier. For the Gibbs test this is based
#' on the ratio of sample range divided by sample standard deviation (Gibbs U-Statistic).
#' }
#' 

re <- "Surv( days_to_death, vital_status )"
va <- setdiff( lst_vars, c( lst_meta, "days_to_last_followup", "days_to_last_known_alive"  )  )

se <- c("gender","years_to_birth")
mod <- fn_train_cox( dat, re, se )
ev  <- fn_eval_cox( dat, re, se, mod, u = 365 )

vck <- check_variables( dat,re, va, min_cases = 35 )
vck <- vck %>%
  mutate( is_cnv = grepl( "_cnv$", variable ), is_cna = grepl( "_cna$", variable ) ) %>% 
  mutate_at(c("check_missing","type","check_entropy"), as.factor )
vck %>% summary

vck %>% 
  ggplot( aes(x=nmiss_pct, y=..density..) ) +
  geom_histogram( bins = 100 ) +
  geom_density( bw = "ucv" ) +
  theme_classic()

vck %>% 
  ggplot( aes(x=entropy, y=..density..) ) +
  geom_histogram( bins = 100 ) +
  geom_density( bw = "ucv" ) +
  theme_classic()

vck %>% 
  ggplot( aes(x=entropy) ) +
  stat_ecdf() + 
  theme_classic()
