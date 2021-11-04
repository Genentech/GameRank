rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )
library( forcats )
library( labelled )
library( tibble )
library( survival )
library( pec )
library( rms )
library( numDeriv )
# library( RccpAlgos )

library( SummarizedExperiment )
library( MultiAssayExperiment)

library( ggplot2 )
library( gridExtra )
library( broom )
library( survival )
library( survminer )

devtools::load_all("~/GameRank")

out_rmd <- "~/GameRank/out"
load( file = "~/GameRank/out/20211103_gmr_result_cox.Rdata" )

selection <- gmr_sel$variable_ranking$variable[1:6]
dt <- bind_rows(
  df_train %>% mutate( ds = 1 ),
  df_eval %>% mutate( ds = 2 ) ) %>%
  bind_rows( df_test %>% mutate( ds = 3 ) )

rmarkdown::render( input = "./R/rmd_game_rank_summary.Rmd", 
                   output_dir = out_rmd,
                   params = list( gmr = gmr_sel ), envir = new.env() )
rmarkdown::render( input = "./R/rmd_calibration_cox.Rmd", 
                   output_dir = out_rmd,
                   params = list( ds = matrix( dt %>% pull( ds ), ncol=1 ), dat = dt, resp = resp1, selection = selection, u = 365, k=1  ), envir = new.env()  )

debugonce(gplot_predictions_cox)
df_train <- params$dat[which(1==params$ds[,params$k]),]
df_eval  <- params$dat[which(2==params$ds[,params$k]),]
df_test  <- params$dat[which(3==params$ds[,params$k]),]
mod <- fn_train_cox( df_train, params$resp, params$selection, u=params$u )
gplot_predictions_cox( df_train, params$resp1, params$selection, mod, u = params$u )
