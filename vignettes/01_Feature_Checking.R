#' ---
#' title: "Feature Checking (Vignette)"
#' output: rmarkdown::html_vignette
#' keep_md: true
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */
# /* knitr::spin(hair="./vignettes/01_Feature_Checking.R", knit=TRUE, report = TRUE, format = "Rmd", precious = TRUE ) */


#' # Introduction
#' 
#' The GameRank package for R intends to provide facilities for building predictive models choosing the best
#' variables for that task. It includes support for regression, logistic regression and survival regression
#' cases by providing standard functions for training and evaluating model performance. As the wrapper approach 
#' allows for direct optimization of measures for the generalization error, a number of combinatorial optimization
#' algorithms for the task of variable selection are included. However wrapper approaches are computationally slow and
#' thus functions for pre-filtering variables to support large-scale model building studies as well as
#' feature constructions functions to define improved non-linear variables and deal with multi-modal distributions 
#' were added.
#' 
#' As Toy Data, these vignettes will make use of four tumor datasets downloaded and prepared from 
#' The Cancer Genome Atlas (TCGA) to exemplify the individual tasks for predictive modelling. Specifically
#' these datasets are for Colon Adenocarcinoma, Diffuse Large B-Cell Lymphoma, Lung Adenocarcinoma and 
#' Uterine Corpus Endometrial Carcinoma. Each dataset includes features taken generated from the CNA and CNV
#' RaggedExperiments through the simplifyTCGA method leading to selection tasks of ~40K features each.
#' 
#' The provided set of vignettes leads through the tasks indicated above but makes use of pre-computed 
#' intermediary outputs to avoid readers waiting for the computational intensive method calls. 
#' 
#' # Building predictive models
#' 
#' In general two distinct settings for predictive modeling can distinguished: 
#' 
#' 1. there are _too many features_ available which leads to a _variable selection_ task, or
#' 2. there are _too few features_ available which leads to a _variable construction_ task.
#' 
#' While the first case can be reasonably approached by algorithmic strategies, the latter leads to an
#' problem with only a limited number of algorithms available and thus requiring mostly creativity and
#' domain knowledge.
#' 
#' For these two settings, there are four steps that can be taken to build a predictive model:
#' 
#' 1. Variable Screening
#' 2. Variable Construction, if necessary
#' 3. Variable Selection, the main task
#' 4. Model validation and checking
#' 
#' Sometimes a fifth step is added to evaluate the robustness of or improve the model, which is
#' 
#' 5. Model exploitation (not discussed here).
#' 
#' # 1. Variable Screening 
#' 
#' As indicated above, an immediate first step is the inspection of available data with respect to:
#' 
#' - _Missing values:_ It's key to understand the level of *missing values* per variable. Multi-variate
#' models may fail to train as the number of complete cases may drop. Alternatively, approaches for
#' handling missing data may be considered. If key variables bear too many missing values, this may
#' be indicative of missing information on a relevant data generating process or even of data quality issues.
#' - _Data type and Information Content:_ Regardless of data type, a variable should comprise a certain level
#' of information, that is should attain more than one value or not be constant except for a few cases. This
#' property can be estimated as the variable *Entropy* and the *Mutual Information* with respect to a response.
#' - _Outliers and Extreme Cases:_ For numerical variables, it is essential to understand the probability and
#' extend of potential *outliers*. Some modeling approaches may be severely impacted by those. We can use two
#' approaches: 1) is based on the so-called _Robust Outlier Test_ using the range outsider [Q1 - c x IQR, Q3 + c x IQR]
#' where c = 1.5 or 3.0 to define outliers, or 2) using the ratio of sample range divided by sample standard deviation, 
#' being the Gibbs U-Statistic. _Remark: The outlier R package provides methods to test for such cases but only for
#' small sample sizes (n<=35) and thus is not used here._
#' 
#' # Using the GameRank package for Variable Screening
#' 
#' Let us start with using GameRank for Variable Screening. First let's load some basic packages and load the data.

#+ eval=TRUE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varchk <- "~/GameRank/out/vig1/"
file_varchk <- file.path( dir_varchk, "tcga_ucec_cna_cnv_varchk.rds" )
file_varchk <- file.path( "~/GameRank/data", "tcga_ucec_cna_cnv_varchk.rds" )

# load data
load( "~/GameRank/data/tcga_ucec_cna_cnv.Rdata" )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )
# Define the outcome and list of variables
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars

#' The code following runs for about 50mins and hence we skip it here and fast forward by loading a results dataset.
#+ eval=FALSE
# Run variable checks
start_time <- Sys.time()
vck <- check_variables( dat, re, va )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Variable checking ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
# Save results and render report
saveRDS( vck, file = file_varchk )
render_variable_checks_summary( vck, dir_varchk )

#+ eval=TRUE
vck <- readRDS( file = file_varchk )
# Summary
vck %>% summary

#' A first glance into the distribution of missing data shows that most variables are set and have 100% non-missing observations.
#' By GameRank package, variables are classified as 'Drop' if the have between 100\% to 30\% missing observations. Variables with 30\%
#' to 20\% are classified as 'Bad', with 20\% to 10\% as 'Try', between 10\% to 1\% as 'Good', and variables with 0\% missing are 
#' classified as 'Perfect'.
#' 
#+ fig.width=7
# Distribution of Missing Data
cut_points <- tribble(
  ~xcut_start, ~xcut_stop, ~xcut_cat,
  0.0,  0.7,  "Drop",
  0.7,  0.8,  "Bad",
  0.8,  0.9,  "Try",
  0.9,  0.99, "Good",
  0.99, 1.0,  "Perfect"
) %>% mutate(xcut_start = xcut_start * 100, xcut_stop = xcut_stop * 100 )
ggplot(  ) +
  geom_rect( aes( xmin = xcut_start, xmax = xcut_stop, ymin=0.0, ymax = 1.0, fill = xcut_cat, color = NA ), data = cut_points, alpha = 0.25 ) +
  geom_histogram( aes(x=p, y=..density..), data = vck, bins = 100 ) +
  geom_density( bw = "ucv" ) + 
  scale_color_brewer( palette = "Accent" ) + 
  xlab( "% non-missing observations")


#' In another evaluation, we can plot the distribution of information content across all variables, measured in bits by the variable entropy.
#' We can see that in this case this distribution is multi-modal with a lower peak around 2.0, another one around 2.3, and two high peaks
#' around 2.65 and 2.75 (orange lines).
#' 
#' Variables with a low entropy do not bear much information and may be dropped in favor of more inforformative ones.
#' 
#+ fig.width=7
# Distribution of Variable Entropies
vck %>% 
  ggplot( aes(x=entropy, y=..density..) ) +
  geom_histogram( bins = 100 ) +
  geom_density( bw = "ucv" ) + 
  geom_vline( xintercept = c(2.0, 2.3, 2.65, 2.75) , color = "orange" )
  
#' Another useful evaluation is the information that a variable shares with regards to the anticipated response variable.
#' This can be expressed as the mutual information, ie the number of bits known about the other variable if only one is
#' known.
#' 
#' While entropy can be used to identify low-information variables, ie those that are constant in the worst case [entropy = 0], the
#' Mutual information can be seen as a generalization of correlation measures. We choose for the information theoretic approach
#' as it can be generalized to fit for continuous, binary and survival variables. In the survival case, the mututal information is
#' estimated with regards to the event flag (as a binary variable).
#' 
#' Entropy and mutual information are estimated on discretized versions of the variable if they are continuous. To this end, the binwidth
#' is determined by leave-one-out cross-validation to obtain an correctly smoothed histogram density estimator. With this binning, the mutual
#' information is then estimated from the counts matrix by the _entropy_ R package
#' 
#' The figure below shows that most included variables do not share more than 4/100 bits about the event rate. This is rather low
#' and might be subject to improvement by variable construction methods in the next vignette.
#' 
#+ fig.width=7
# Distribution of Mutual Information with Response
vck %>% 
  ggplot( aes(x=mutual_information, y=..density..) ) +
  geom_histogram( bins = 100 ) +
  geom_density( bw = "ucv" )

#' Putting entropy and mutual information together provides a comprehensive view about variability and information with respect
#' to the outcome. We can see here that even high entropy variables do not necessarily share most information about the desired
#' outcome, while the mid-range entropy variables share most.
#' 
#' This plot suggesting a selection of variables above 0.04 for mutual information plus for exploration some that are >=0.02 and 
#' have entropy >3 if they do not bear missing data..
#' 
#+ fig.width=7
# Joint distribution of variable entropy and mutual information with response
vck %>% 
  ggplot( aes(x=entropy, y=mutual_information, color = check_missing ) ) +
  geom_point()


#' A final consideration should be given to outliers and thus to variable skew. Outliers can be on both
#' ends of a continuous scale. We can see that many of the informative points do have outliers, which may
#' distort the model. Nontheless we should be exploring them but also see if there are variables with
#' not too many outliers on each end that can be used.
#' 
#+ fig.width=7
# Scatter Plot of Robust Outliers Min and max colored by Entropy
vck %>%
  mutate( flag = (mutual_information > 0.04 | (mutual_information > 0.02 & entropy > 3)) ) %>%
  ggplot( aes( x=rot.nmin, y=rot.nmax, color=flag ) ) +
  geom_point() +
  scale_color_discrete()


#' If we factor the mutual information with respect to event rate in, we can see that variables with a low
#' number of outliers at the lower end have better mutual information than vice versa. 
#' 
#' Thus, we let's keep also variables that have <10 outliers either at each end.
#' 
#+ fig.width=7
# Scatter Plot of Robust Outliers Min and max colored by Mutual Information
vck %>%
  ggplot( aes( x=rot.nmin, y=rot.nmax, color=mutual_information ) ) +
  geom_point()

#' Our final list of candidate variables thus might be, after adjusting the mutual information to 0.045, and
#' entropy to 3.1, is giving us 279 variables for model building via wrapper algorithms:
ret <- vck %>% 
  filter( "Perfect"==check_missing ) %>% 
  filter( (mutual_information > 0.045 | (mutual_information > 0.02 & entropy > 3.1) | (rot.nmin < 10 & rot.nmax < 10) ) )
vars <- ret %>% pull( variable ) %>% unique
cat( "List of choosen variables : \n")
vars

#' Without domain knowledge, this selection process is more or less arbitrary. It definitely should be repeated and also include
#' an random element such that the set of variables is not too much biased towards the one or the other feature and also allows for
#' completely odd features to be selected. Always keep in mind that at the end, the *combination* of variables may lead to a successful model
#' and not its individual components alone.
