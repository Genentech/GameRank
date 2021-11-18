#' ---
#' title: "GameRank Algorithm (Vignette)"
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

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varchk <- "out/vig2/"
file_varchk <- file.path( dir_varchk, "tcga_dlbcl_cna_cnv_construction.rds" )

#' Example code
load("data/tcga_dlbcl_cna_cnv.Rdata")
lst_tumor_genes <- readRDS("data/cancer_genes.rds")

# Construct new features for selection
# they are just string expressions
ll <- setdiff( lst_vars, lst_clinvars )
tib <- dplyr::inner_join(
  tibble( var_cna = grep( "_cna$", ll, value = TRUE ) ) %>% mutate( var = gsub( "_cna$", "", var_cna ) ),
  tibble( var_cnv = grep( "_cnv$", ll, value = TRUE ) ) %>% mutate( var = gsub( "_cnv$", "", var_cnv ) ),
  "var" )
tib <- tib %>% 
  mutate( feat_ratio = sprintf( "( %s / %s )", var_cna, var_cnv ),
          feat_atan  = sprintf( "atan( %s / %s )", var_cna, var_cnv ) )
# Let's annotate the genes as if they are know as tumor suppressors or proto-oncogenes
tib <- tib %>% mutate( is_tumor_suppressor = (var %in% lst_tumor_genes$tumor_suppressors),
                       is_proto_oncogene   = (var %in% lst_tumor_genes$proto_oncogenes) ) %>%
  mutate( is_cancer_gene = is_tumor_suppressor | is_proto_oncogene )
lst_constr <- union( tib %>% pull( feat_ratio ), tib %>% pull( feat_atan ) )

# There is a set of 6 genes that are both: tumor suppressor and proto-oncogenes per UniProtKB as of 2021-11-17
tibx <- tib %>% filter( is_tumor_suppressor & is_proto_oncogene ) %>% as.data.frame 
datx <- dat %>% dplyr::select( all_of( c(lst_keys, lst_meta, lst_outcomes, tibx %>% pull( var_cna ), tibx %>% pull( var_cnv ) )) )
datx

# TODO Construct Box-Cox for age and some tumor suppressors
datx %>%
  ggplot( aes( x=(MAF_cna / MAF_cnv)^(+3.0), y=..density.. ) ) +
  geom_histogram( bins = 10 ) +
  geom_density( bw = "ucv" )

# debugonce(simple_transforms)
trdat <- simple_transforms( datx, tibx %>% pull(var_cna) )
trdat$transformations %>% Reduce( dplyr::bind_rows, ., NULL ) %>% group_by( variable ) %>% filter( max(W)==W )

# There is no Box-Cox transform yet for Survival models, hence we demonstrate
# it for predicting the events component through age
# debugonce(box_cox_binomial)
bcdat <- box_cox_binomial( dat %>% mutate( vital_status = as.logical( vital_status ) ) , 
                           "vital_status", c("years_to_birth", tibx %>% pull( feat_ratio )) )
bcdat$transforms %>% map_dfr( function( dd ) { dd[["boxcox_result"]] <- NULL; return( as_tibble( dd ) )  } ) %>% as.data.frame

bcdat$data %>%
  ggplot( aes(x=years_to_birth_bc_bin, y=..density.., group=vital_status, color=vital_status ) ) +
  geom_histogram() + geom_density( bw = "ucv" )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )
re <- "Surv( days_to_death, vital_status )"

# Run variable checks
start_time <- Sys.time()
vck <- check_variables( dat, re, va )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Variable checking ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )

# Save results and render report
saveRDS( vck, file = file_varchk )
render_variable_checks_summary( vck, dir_varchk )
