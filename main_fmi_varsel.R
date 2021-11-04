rm( list=ls() )
# Variable selection ----
# Meta variables
file_var_rds <- "~/GameRank/tmp/df_prima_selection_data.rds"
lst_keys <- c("STUDYID",
              "USUBJID",
              "A_RND", #  Randomized Patients
              "A_MEXT", # Maintenance Extent Of Exposure (DAYS)
              "A_MAX" # Maximum Observation Time 
)

# Response variables
lst_resp <- c(
  "OS_AVAL","OS_EVENT", 
  "PFSINV_AVAL", "PFSINV_EVENT",
  "A_IAL", # Induction Patient Alive
  "A_IDD", # Induction Patient Death
  "A_MAL", # Maintenance Patient Alive
  "A_MDD"  # Maintenance Patient Death
)

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

devtools::load_all("~/GameRank")

df_sel <- readRDS( file = file_var_rds )
df_sel %>% summary

lst_vars <- df_sel %>% colnames %>% setdiff( c(lst_keys, lst_resp) )
lst_is_char <- sapply( lst_vars, function( v ) is.character(df_sel[[v]]) )
df_sel <- df_sel %>% 
  mutate_at( names(which(lst_is_char)), function(x) x %>% as.factor %>% fct_drop ) %>%
  mutate( A_IDD = ("YES"==A_IDD) )
df_sel %>% summary

g_set <- rep_len( c(1,1,1,2,3), length.out = nrow(df_sel) )
g_set <- g_set[order( runif( length(g_set) ) )]
table(g_set)
prop.table(table(g_set))

var_chks <- check_variables( df_sel, NULL, lst_vars )
lst_vars <- var_chks %>% arrange( desc( entropy )  ) %>% pull( variable ) # Sort variables by entropy to ensure worst informative variables are removed first


df_sel <- df_sel %>% mutate_at( names(which(lst_is_char)), function(x) x %>% as.factor %>% fct_drop )
df_test  <- df_sel[which(3==g_set),]
df_eval  <- df_sel[which(2==g_set),]
df_train <- df_sel[which(1==g_set),]
df_trnevl <- bind_rows( df_train, df_eval )

ds_prior <- prepare_splits( ds = 3L, dat=df_trnevl, resp = resp1, vars = lst_vars, fn_train = fn_train_cox, fn_eval = fn_eval_cox ) %>% as.matrix


resp1 <- "Surv( PFSINV_AVAL, PFSINV_EVENT )"
# resp1 <- "A_IDD"
gmr_sel <- game_rank( dat = df_trnevl, resp = resp1, vars = lst_vars, fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 7L, team_size = 20L, rounds = 25L, maximize = FALSE, u = 365 )
gmr_sel
save( file_var_rds, lst_keys, lst_resp, gmr_sel, df_sel, df_test, df_train, df_eval, df_trnevl, lst_vars, resp1, file = "~/GameRank/out/20211103_gmr_result_cox.Rdata" )


# rnd_sel <- random_selection( dat = df_trnevl, resp = resp1, vars = lst_vars, fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 7L, ds = ds_prior, nevals = 100L, maximize = FALSE, u = 365 )
# rnd_sel
# 
# fwd_sel <- forward( dat = df_trnevl, resp = resp1, vars = lst_vars, fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 2L, ds = ds_prior, maximize = FALSE, u = 365 )
# fwd_sel

# bwd_sel <- backward( dat = df_trnevl, resp = resp1, vars = lst_vars[1:5], fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 3, ds = ds_prior, maximize = FALSE, u = 365 )
# bwd_sel

# bwd_sel <- bidirectional( dat = df_trnevl, resp = resp1, vars = lst_vars[1:10], fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 3, ds = ds_prior, maximize = FALSE, u = 365 )
# bwd_sel

# lrs_sel <- lrsearch( dat = df_trnevl, resp = resp1, vars = lst_vars[1:10], fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 3L, L = 3L, R = 2L, ds = ds_prior, maximize = FALSE, u = 365 )
# lrs_sel

#
# To Do List:
#  - Add Oscillating Search Algorithm for Feature Selection
#  - Debug floating searches
# TODO Debug and fix
# sffs_sel <- sffs( dat = df_trnevl, resp = resp1, vars = lst_vars[1:10], fn_train = fn_train_cox, fn_eval = fn_eval_cox, m = 3L, ds = ds_prior, maximize = FALSE, u = 365 )
# sffs_sel
#

