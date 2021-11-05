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

library( numDeriv )
library( pec )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_coad_cna_cnv.Rdata")
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
va <- setdiff( lst_vars, c( lst_meta, "days_to_last_followup" )  )

se <- c("gender","years_to_birth")
mod <- fn_train_cox( dat, re, se )
ev  <- fn_eval_cox( dat, re, se, mod, u = 365 )

vck <- check_variables( dat,re, va )
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

var_cnv <- vck %>% filter( 100.0==p ) %>% filter( is_cnv ) %>% filter( max(entropy) * 0.975 <= entropy ) %>% pull( variable )
var_cna <- vck %>% filter( 100.0==p ) %>% filter( is_cna ) %>% filter( max(entropy) * 0.75 <= entropy ) %>% pull( variable )
intersect( gsub( "_cna$", "", var_cna ), gsub( "_cnv$", "", var_cnv ) )
var_clin <- vck %>% filter( !is_cnv & !is_cna ) %>% filter( check_missing %in% c("Perfect","Good")) %>% pull( variable )
sel_vars <- union( var_clin, var_cnv ) %>% union( var_cna )

ds <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) )
ds <- ds[ order( runif( length(ds) ) ) ]

df_train <- dat[which(1==ds),]
df_eval  <- dat[which(2==ds),]
df_test  <- dat[which(3==ds),]
df_trev <- bind_rows( df_train, df_eval )

grm <- game_rank( df_trev, re, sel_vars, fn_train_cox, fn_eval_cox, m = 5, maximize = FALSE, team_size = 30L, rounds = 10L, min_matches_per_var = 5L, u = 365 )
save.image( file = "~/GameRank/out/20211105_gmr_coad.rds")

