#' ---
#' title: "Feature Construction (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */


#' 
#' In the last vignette, we realized that the raw CNA and CNV features do not bear more than 4/100 bits of information about the overall survival rate.
#' This seems rather low and certainly warrants the attempt to improve by feature construction.
#' Within this vignette we'll outline how this can be done:
#' 
#' - First, we'll try by trying a functional relationship, that is evaluating the ratio of CNA to CNV for a given gene as well as the arc tanget of CNA and CNV.
#' The motivation for the arc tangent is that it never attains invalid values (using the atan2 R function) and yields values within +1.0 and -1.0 which allow model coefficients
#' to stay within a reasonable range.
#' - Secondly, we'll demonstrate applying the 'usual transformations' - log, sqrt, cube root and z-score - to achieve better normality on a small subset of features stemming from 
#' combined tumor suppressor and proto-oncogenic genes.
#' - Finally, we'll exercise using the Box-Cox transformation to obtain power transformed variables with respect to binomial outcomes (in this case the event rate).
#' 


#' Initialization and setup code
#+ eval=TRUE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varchk <- "~/GameRank/out/vig2/"
file_varchk <- file.path( dir_varchk, "tcga_ucec_cna_cnv_construction.rds" )
file_varchk <- file.path( "~/GameRank/data", "tcga_ucec_cna_cnv_construction.rds" )
file_smtf <- file.path( dir_varchk, "tcga_ucec_cna_cnv_simple_transforms.rds")
file_smtf <- file.path( "~/GameRank/data", "tcga_ucec_cna_cnv_simple_transforms.rds" )

load( "~/GameRank/data/tcga_ucec_cna_cnv.Rdata" )
lst_tumor_genes <- readRDS("~/GameRank/data/cancer_genes.rds")

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )
re <- "Surv( days_to_death, vital_status )"

#' # Construct new variables and check them
#' Next let's make a list of potential features that may improve. As features are handled as character strings
#' in GameRank and fed via the formula interface into the training:validation functions, they can be written as
#' formula expressions, including - like for atan2 - function calls.
#'
#' Thus we merge the variables for CNA and CNV by gene name, and construct strings for computing the ratio
#' and the atan2( CNA, CNV ) functions as features. Those will be screened next by the variable checks procedure.
ll <- setdiff( lst_vars, lst_clinvars )
tib <- dplyr::inner_join(
  tibble( var_cna = grep( "_cna$", ll, value = TRUE ) ) %>% mutate( var = gsub( "_cna$", "", var_cna ) ),
  tibble( var_cnv = grep( "_cnv$", ll, value = TRUE ) ) %>% mutate( var = gsub( "_cnv$", "", var_cnv ) ),
  "var" )
tib <- tib %>% 
  mutate( feat_ratio = sprintf( "( %s / %s )", var_cna, var_cnv ),
          feat_atan  = sprintf( "atan2( %s , %s )", var_cna, var_cnv ) )
# Let's annotate the genes as if they are know as tumor suppressors or proto-oncogenes
tib <- tib %>% mutate( is_tumor_suppressor = (var %in% lst_tumor_genes$tumor_suppressors),
                       is_proto_oncogene   = (var %in% lst_tumor_genes$proto_oncogenes) ) %>%
  mutate( is_cancer_gene = is_tumor_suppressor | is_proto_oncogene )
lst_constr <- union( tib %>% pull( feat_ratio ), tib %>% pull( feat_atan ) )
va <- lst_constr

#' The following code runs for about 50 mins.
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

#+ fig.width=7
vck %>% 
  ggplot( aes(x=entropy, y=mutual_information, color = check_missing ) ) +
  geom_point()

#' # Evaluate the standard transformations (square root, cube root, log and z-score)
#' In the first attempt, there was no observable increase in mutual information seen. While there exist two clusters of which one shows an increase in
#' entropy around 3.0, none of the transformed variables gave an mutual information exceeding 0.04. Based on these, results continuing with the originally 
#' selected set of features seems the right thing to do.
#' 
#' However, there are still options to pursue: we can apply feature construction methods or strategies. The basic one is to try a set of standard transforms
#' and check if they improve normality of the distribution. Normality can be assessed by the W-statistic of the Shapiro-Wilk test. An increase in the W-statistic
#' indicates improved Normality.
#' 
#+ eval=FALSE
start_time <- Sys.time()
datsmp <- simple_transforms( dat, lst_vars )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Simple transformations ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
tbsmp <- datsmp$transformations %>% Reduce( bind_rows, ., NULL ) %>% as_tibble
saveRDS( tbsmp, file = file_smtf )

#+ eval=TRUE
tbsmp <- readRDS( file = file_smtf )

#' We have checked our original list of variables now for whether the usual transformation (square root, cube root, log, z-score) improve their Normality.
#' Let's have a look how well which transformation performs:
tbbst <- tbsmp %>% group_by( variable ) %>% filter( max(W)==W )
tbbst$transform %>% table 
tbbst$transform %>% table %>% prop.table %>% sprintf( "%1.4f", . )

#' Transformations in nearly all cases improve normality. In 53.4\% a cube root transform
#' followed by a z-score transform in 43\% works best. The log transform helps in 
#' about 3\% and the square root in just 0.01\%. The raw variable is best just in 0.09\%
#' of the cases.
#' 
#' 

#'
#' # Construct power-transformed variables by the Box-Cox transformation
#'  
#' In some cases, non-linear features can be constructed using the Box-Cox transform. This approach estimates a maximum-likelihood parameter p for the univariate model Y ~ x^p.
#' Currently, this approach is only available for the regression and binomial/logistic regression case.
#' 
#' We'll outline it on a small sample of selected variables, namely those for which the proteins were flagged as both: tumor suppressor and proto-oncogenes in UniProtKB.
#'
#' There is a set of 5 genes that are both: tumor suppressor and proto-oncogenes (per UniProtKB as of 2021-11-17) for which we'll evaluate a Box-Cox transformation
#+ eval=TRUE
tibx <- tib %>% filter( is_tumor_suppressor & is_proto_oncogene ) %>% as.data.frame 
tibx
datx <- dat %>% dplyr::select( all_of( c(lst_keys, lst_meta, lst_outcomes, lst_clinvars, tibx %>% pull( var_cna ), tibx %>% pull( var_cnv ) )) )
datx

#' As there is no Box-Cox transformation algorithm defined yet for Survival models, we'll demonstrate
#' it for predicting the events component through the set of 6 CNA/CNV variables
bcdat <- box_cox_binomial( dat %>% mutate( vital_status = as.logical( vital_status ) ) , 
                           "vital_status", c(tibx %>% pull( var_cna ), tibx %>% pull( var_cnv ) ) )
dfbcres <- bcdat$transforms %>% map_dfr( function( dd ) { dd[["boxcox_result"]] <- NULL; return( as_tibble( dd ) )  } ) %>% as.data.frame
dfbcres

#' In our variable selection process, we can now integrate the terms provided. Columns py denotes the odds ratio and px denotes the exponent for the
#' univariate model predicting vital status.
#' 


#' # Check for multi-modal variables
#' 
#' Since there is no age or other known multi-modal feature available, let's pretend that the days to death are a feature and show how
#' the search for multi-modal distributions can be done in GameRank.
#' 
#' Evaluation of multi-modal distributions is done by fitting a Gaussian Mixture Model with up to n_comp components. For each number of components
#' m_fits models are fit (using flexmix) to estimate the best Akaike Information Criterion. The final number of components is determined then as
#' the minimal AIC with the minimal number of components for all models for which at least min_prop_converged models were obtained.
#' 
mumo <- check_multimodality( dat, NULL, c("days_to_death"), n_comp = 5, m_fits = 5  )
mumo$transforms
print( mumo$transforms$days_to_death$best_model )
pm <- flexmix::parameters( mumo$transforms$days_to_death$best_model )
pm
po <- flexmix::prior( mumo$transforms$days_to_death$best_model ) 
po
mumo$transforms$days_to_death$cut_points

#' In this case we obtain a mixture distribution of two components: one for the censored patients around 3424 days and one for the patients
#' dying during the early days of observation. Check for multi-modality provides a numeric array for cut-points. These were found by identifying
#' the intersection of the component distributions multiplied by their prior.
#' 
#' We can plot the mixture distribution then:
#+ eval=FALSE, fig.width=7
# TODO Find better variable as showcase
pt <- tibble( x = dat$days_to_death )
di <- tibble( x = seq( min(dat$days_to_death), max(dat$days_to_death), 0.1 ) ) %>%
  mutate( di1 = dnorm( x, mean = pm[1,1], sd = pm[2,1] ) * po[1],  # Calculate distribution for component 1 with prior 1
          di2 = dnorm( x, mean = pm[1,2], sd = pm[2,2] ) * po[2] ) # Calculate distribution for component 2 with prior 2
ggplot() +
  geom_histogram( aes(x=x,y=..density..), data = pt, bins = 100 ) +
  geom_density(  aes(x=x,y=..density..), data = pt, bw = "ucv", color = "orange" ) +
  geom_line( aes( x=x, y=di1 ), data = di, color = "red" ) +
  geom_line( aes( x=x, y=di2 ), data = di, color = "green" )
