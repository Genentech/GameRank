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

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

#' Example code
load("data/tcga_dlbcl_cna_cnv.Rdata")
dat <- dat %>% mutate( days_to_death = ifelse( !is.na(days_to_death), days_to_death, 1+max(days_to_death, na.rm=TRUE)))

#' # Introduction
#' 
#' This Vignette aims to introduce feature selection using the GameRank algorithm.
#' We'll use the Colon Adenocarcinoma data from TCGA with features extracted for
#' Copy Number Aberrations and Copy Number Variations. This dataset has N=457 patients
#' included and our goal will be to find a model predicting Overall Survival.
#+ fig.width=7
fit <- survfit( Surv( days_to_death, vital_status ) ~ 1, data = dat )
ggsurvplot(fit,dat)

re <- "Surv( days_to_death, vital_status )"
va <- setdiff( lst_vars, c( lst_meta, "days_to_last_followup","colname" )  )
va <- va[1:25]

se <- c("gender","years_to_birth")
mod <- fn_train_cox( dat, re, se )
ev  <- fn_eval_cox( dat, re, se, mod, u = 365 )

fn_bt_eval_cox <- ff_fn_eval_bootstrap( fn_train_cox, fn_eval_cox, 10L )
mod <- fn_train_dummy( dat, re, se )
bev  <- fn_bt_eval_cox( dat, re, se, mod, u = 365 )

# mod <- fn_train_dummy(dat,re,se)
# debugonce(fn_bt_eval_cox)
# bev  <- fn_bt_eval_cox( dat, re, se, mod, u = 365 )

# ds <- GameRank::build_splits(2, dat, re, se, fn_train_dummy, fn_bt_eval_cox, u = 365 )
# fsl <- forward( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 5, ds, FALSE, u = 365)
# saveRDS( fsl, file = "~/tmp/20211108_forward.rds" )

# debugonce(fn_bt_eval_cox)
gmr <- game_rank( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 5L, c(2L,2L), FALSE, 3L, 10L, 5L, u = 365  )
gmr


fsl <- game_rank( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 5, FALSE, u = 365)
fsl


bsl <- backward( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 5, ds, FALSE, u = 365)
bsl
