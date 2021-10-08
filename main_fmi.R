rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )
library( labelled )

library( SummarizedExperiment )
library( MultiAssayExperiment)

devtools::load_all("~/GameRank")

src.wd <-  "/opt/bee/analyses/CIN/cin_1269/lymphoma/prima/"

data.mae <- readRDS( file.path( src.wd,  "mae.rds" ) )
data.mae

# Build dataset ----
# Meta variables
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
# Predictor variables
lst_vars <- c(
  "A_IREG", #   Induction Regimen
  "A_IRSP", #   Response To Induction Regimen
  "A_CNTRY", # Nation Or Country
  
  "RND1", #  Randomized treatment (Maintenance phase)
  "AGE",
  "SEX",
  "WGTKG",
  "RACE",
  "TOBAYN",
  "ALCYN",
  "CAFYN",
  # "REPSTAT", #  Female reproductive status
  "BDSURF",
  
  "A_HGBC", #   Categorized Hemoglobin (G/DL)
  "A_HGBN", #   Categorized Hemoglobin, Num
  "A_HGBCM", #   BS.M. Categorized Hemoglobin (G/DL)
  "A_HGBNM", #   BS.M. Categorized Hemoglobin, Num
  "A_LDHC", #   Categorized LDH (UI/L)
  "A_LDHN", #   Categorized LDH, Num
  "A_LDHCM", #   BS.M. Categorized LDH (UI/L)
  "A_LDHNM", #   BS.M. Categorized LDH, Num
  "A_BT2C", #   Categorized Beta2 Microglobulin (MG/L)
  "A_BT2N", #   Categorized Beta2 Microglobulin, Num
  "A_BT2CM", #   BS.M. Categorized Beta2 Microglob.(MG/L)
  "A_BT2NM", #   BS.M. Categorized Beta2 Microglob., Num
  "A_ANN", #   Ann Arbor Stage
  "A_BINV", #   Bone Marrow Involvement
  "A_BINVN", #   Bone Marrow Involvement, Num
  "A_BINVM", #   BS.M. Bone Marrow Involvement
  "A_BINVNM", #   BS.M. Bone Marrow Involvement, Num
  "A_BSYM", #   B Symptoms
  "A_BSYMN", #   B Symptoms, Num
  "A_BSYMM", #   BS.M. B Symptoms
  "A_BSYMNM", #   BS.M. B Symptoms, Num
  "A_ECOG", #   ECOG
  "A_ECOGM", #   BS.M. ECOG
  "A_EXTC", #   No. Of Reported Extra-Nodal Sites
  "A_EXTA", #   Extra-Nodal Involvement 
  "A_EXTD", #   No. Of Derived Extra-Nodal Sites
  "A_PERC", #   No. Of Reported Peripheral Lymph Sites
  "A_PERD", #   No. Of Derived Peripheral Lymph Sites
  "A_AREA", #   No. Of Nodal Areas
  
  "A_FCRF", #  FLIPI (Reported On CRF Page 7)
  "A_FDER", #  FLIPI (Derived From CRF Pages)
  
  "A_PROG", #  Baseline Prognosis
  # "A_PROGN", #  Baseline Prognosis, Num (?)
  
  "A_HTBC", #  High Tumor Burden GELF Criteria
  "A_HTBN", #  High Tumor Burden GELF Criteria, Num
  
  "A_BULDC", #  Bulky Disease
  "A_BULDN", #  Bulky Disease, Num
  "A_GEBSC", #  B Symptoms
  "A_GEBSN", #  B Symptoms, Num
  
  "A_GLBC", #  Elev. Serum LDH Or Beta2Microglob.
  "A_GLBN", #  Elev. Serum LDH Or Beta2Microglob., Num
  "A_INSC", #  Involv. Of At Least 3 Nodal Sites
  "A_INSN", #  Involv. Of At Least 3 Nodal Sites, Num
  "A_SPEC", #  Splenic Enlargement
  "A_SPEN", #  Splenic Enlargement, Num
  "A_CSYC", #  Compressive Syndrome
  "A_CSYN", #  Compressive Syndrome, Num
  "A_PLPC", #  Pleural Peritoneal Effusion
  "A_PLPN", #  Pleural Peritoneal Effusion, Num
  
  "A_IPD", #  Initial Pathology Diagnosis
  "A_LPD", #  Local Pathology Diagnosis
  "A_IDG", #  Independent Pathology Diagnosis
  "A_IPCF", #  Initial Pathology Diagnosis Confirm.
  "A_IPCFN"  #  Initial Pathology Diagnosis Conf., Num
  
  # "ALBUM_LBUNIT", #  Albumin lab unit
  # "ALBUM_LBVAL", #  Albumin lab value
  # 
  # "B2MICRO_LBUNIT", #  B2-Microglobulin lab unit
  # "B2MICRO_LBVAL", # B2-Microglobulin lab value
  # 
  # "HGB_LBUNIT", #  Hemoglobin lab unit
  # "HGB_LBVAL", #  Hemoglobin lab value
  # 
  # "LDH_LBUNIT", #  LDH lab unit
  # "LDH_LBVAL", #  LDH lab value
  # 
  # "PLATE_LBUNIT", #  Platelet lab unit
  # "PLATE_LBVAL"   #  Platelet lab value
  )
length(lst_vars)


info <- tibble( var_names = colData(data.mae) %>% colnames, 
                var_info  = colData(data.mae) %>% as_tibble %>% var_label ) %>% unnest(var_info)
# info %>% View


dat <- colData(data.mae) %>% 
  as_tibble %>%
  dplyr::select( all_of( c( lst_keys, lst_resp, lst_vars ) ) )
dat %>% head
dat %>% str
