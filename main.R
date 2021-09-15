rm( list=ls() )
library( dplyr )
library( tidyr )
library( tibble )
library( purrr )

library( survival )
library( rms )
library( pec )
library( RcppAlgos )
library( labelled )

library( yaml )
library( SummarizedExperiment )
library( MultiAssayExperiment )

yaml_file <- "~/CoMMpass_Dashdis_Data/CoMMpass-integration-dashdis.yaml"
yaml_cfg <- read_yaml( file = yaml_file )

library( numDeriv )
library( formula.tools )
library( ROCR )

source( "R/model_functions.R" )
source( "R/game_rank.R" )
source( "R/backward.R" )
source( "R/forward.R" )
source( "R/random.R" )
source( "R/variable_construction.R" )
source( "R/descriptive_summary.R" )

# Load CoMMpass Labs data as test set for variable selection
path.dashdis_data <- "~/CoMMpass_Dashdis_Data/"
file_adsl <- yaml_cfg$output$file$adsl
df_adsl <- readRDS( file = file.path( path.dashdis_data, file_adsl ) )
df_adsl
file_labs <- yaml_cfg$output$file$se$labs
se_labs <- readRDS( file = file.path( path.dashdis_data, file_labs ) )
se_labs
full_data <- as_tibble( colData( se_labs ) )
full_data %>% str

df_lab <- as.data.frame( t( assays(se_labs)[['lab']] ) ) %>%
  tibble::rownames_to_column("SPECTRUM_SEQ") %>%
  as_tibble %>%
  dplyr::left_join( as_tibble(colData(se_labs)), "SPECTRUM_SEQ" ) %>%
  dplyr::select( SPECTRUM_SEQ, USUBJID, tidyr::everything() ) %>%
  dplyr::left_join( df_adsl, "USUBJID" ) 
str(df_lab)
df_lab[1:10,1:10]

df_lab %>% colnames %>% dput
resp_cols <- c(
  "DEATH_AVALC", 
  "PFS_AVAL", "PFS_EVENT", 
  "PFSLN1_AVAL", "PFSLN1_EVENT", 
  "PFSLN2_AVAL", "PFSLN2_EVENT", 
  "PFSLN3_AVAL", "PFSLN3_EVENT", 
  "OS_AVAL", "OS_EVENT"
)
key_cols <- c("SPECTRUM_SEQ", "USUBJID")
sel_cols <- c("SPECTRUM_SEQ", "USUBJID",
              
              "BL_MPROTEIN", 
              "DEATH_AVALC", 
              "PFS_AVAL", "PFS_EVENT", 
              "PFSLN1_AVAL", "PFSLN1_EVENT", 
              "PFSLN2_AVAL", "PFSLN2_EVENT", 
              "PFSLN3_AVAL", "PFSLN3_EVENT", 
              "OS_AVAL", "OS_EVENT", 
              
              "D_LAB_cbc_abs_neut", "D_LAB_chem_albumin", 
              # "D_LAB_chem_bun",
              "D_LAB_chem_calcium", "D_LAB_chem_creatinine", 
              "D_LAB_chem_glucose", "D_LAB_cbc_hemoglobin", "D_LAB_serum_kappa", 
              "D_LAB_chem_ldh", "D_LAB_serum_m_protein", "D_LAB_cbc_platelet", 
              "D_LAB_chem_totprot", "D_LAB_cbc_wbc", "D_LAB_serum_iga", "D_LAB_serum_igg", 
              "D_LAB_serum_igm", "D_LAB_serum_beta2_microglobulin", "D_LAB_serum_lambda", 
              "D_LAB_urine_24hr_m_protein", "D_LAB_urine_24hr_total_protein", 
              "D_LAB_serum_c_reactive_protein", "D_LAB_serum_igd", "D_LAB_serum_ige", 
              
              "D_PT_raceoth", "D_PT_race", "D_PT_ethnic", "D_PT_gender", 
              "IMWG_Risk_Class", 
              "Prolif_Index", 
              "D_PT_age", 
              "demog_height", 
              "demog_weight", "D_PT_iss", 
              "DEMOG_HEIGHTUNITOFM", "DEMOG_WEIGHTUNITOFM", 
              
              "D_PT_issstage_char", "D_PT_issstage", "CREATININE", "ecog", 
              "sct_elig",
              
              # "D_PT_therclass", "D_PT_therfstn", "D_PT_therclassn", 
              # "D_PT_maxline", "ftrttrpl", "sct_bresp", "line1sct", "sctline", 
              # "prescttrt", "presctst", "prescten", "fstsctdy",
              
              # "bmtx_type", 
              # "sctflag", "sct1stf", "mainttrt", "maintstdy", "maintendy", "maintday", 
              # "maintflag", "maint1stf", "D_PT_pddy", "D_PT_pdflag", "D_PT_ttfpdw", 
              # "D_PT_respdur", 
              # "D_PT_dresp", "earlydth", "fhr", 
              
              "D_PT_mmstatus", "D_PT_mmstatus1", "D_PT_mmstatus2", 
              "D_PT_mmstatus3", 
              "D_PT_rapd", 
              
              # "demog_vj_interval", "demog_visitdy",
              
              "DEMOG_PATIENTAGE", "DEMOG_GENDER", 
              "DEMOG_AMERICANINDIA", "DEMOG_ASIAN", "DEMOG_BLACKORAFRICA", 
              "DEMOG_NATIVEHAWAIIA", "DEMOG_WHITE", "DEMOG_OTHER", "DEMOG_SPECIFY", 
              "DEMOG_ETHNICITY"
              
              # "DEMOG_DAYOFVISIT", "demog_visit", "STUDY_ACRONYM"
)
df_lab <- df_lab %>% 
  dplyr::select( all_of(sel_cols) ) %>%
  mutate_if( is.character, as.factor )

# Run forward selection
resp <- "DEATH_AVALC"
resp <- "D_PT_age"
resp <- "Surv(OS_AVAL,OS_EVENT)"
vars <- c("D_LAB_cbc_abs_neut", "D_LAB_chem_albumin", 
          # "D_LAB_chem_bun", 
          "D_LAB_chem_calcium", "D_LAB_chem_creatinine", 
          # "D_LAB_chem_glucose",
          "D_LAB_cbc_hemoglobin", "D_LAB_serum_kappa", 
          # "D_LAB_chem_ldh", 
          "D_LAB_serum_m_protein", "D_LAB_cbc_platelet", 
          "D_LAB_chem_totprot", "D_LAB_cbc_wbc", 
          # "D_LAB_serum_iga",
          "D_LAB_serum_igg", 
          # "D_LAB_serum_igm", 
          # "D_LAB_serum_beta2_microglobulin", 
          "D_LAB_serum_lambda"
          # "D_LAB_urine_24hr_m_protein"
          # "D_LAB_urine_24hr_total_protein"
          # "D_LAB_serum_c_reactive_protein",
          # "D_LAB_serum_igd", 
          # "D_LAB_serum_ige"
          )
summary( df_lab[,vars] )
dd <- df_lab[which(complete.cases(df_lab[,vars])),vars]
summary(dd)
nrow(dd)
dd <- df_lab[which(complete.cases(df_lab[,vars])),]

debugonce(simple_transforms)
rr <- simple_transforms( dat =dd, vars = vars )



# fres <- forward( dat = dd, resp = resp, vars = vars,
#                  fn_train = fn_train_cox, fn_eval = fn_eval_cox, u = 365.25, maximize = TRUE )
# fres
# bres <- backward( dat = dd, resp = resp, vars = vars,
#                   fn_train = fn_train_cox, fn_eval = fn_eval_cox, u = 365.25, maximize = TRUE  )
# bres
# gres <- game_rank( dat = dd, resp = resp, vars = vars,
#                    fn_train = fn_train_cox, fn_eval = fn_eval_cox, u = 365.25, maximize = FALSE  )
# gres
# rres <- random_selection( dat = dd, resp = resp, vars = vars,
#                           fn_train = fn_train_cox, fn_eval = fn_eval_cox,
#                           ksize = 3L, min_sample_per_var = 4L, u = 365.25, maximize = FALSE )
# rres


# vars <- c( "D_LAB_cbc_hemoglobin", "D_LAB_serum_kappa" )
# mo <- fn_train_cox( dat = dd, resp, vars )
# debugonce(fn_eval_normal)
# va <- fn_eval_cox( dat = dd, resp, vars, mo, u = 365.25 )
