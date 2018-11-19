# Config


library(tidyverse)
source('config.R')


sim_score_av <- read_csv(file.path(base_dir,"sim_matrix_resnik_avg_314_log2_2018-11-11.csv")) %>% 
  column_to_rownames('X1')

sim_score_max <- read_csv(file.path(base_dir,"sim_matrix_resnik_def_314_log2_2018-11-09.csv")) %>% 
  column_to_rownames('X1')

n1_100k_av <- read_csv(file.path(base_dir,"100k_2pats_resnik_avg_314_log2_v2_2018-11-12.csv")) 
n3_100k_av <- read_csv(file.path(base_dir,'100k_3pats_resnik_avg_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n4_100k_av <- read_csv(file.path(base_dir,'100k_4pats_resnik_avg_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n5_100k_av <- read_csv(file.path(base_dir,'100k_5pats_resnik_avg_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n7_100k_av <- read_csv(file.path(base_dir,'100k_7pats_resnik_avg_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))

n1_100k_max <- read_csv(file.path(base_dir,"100k_2pats_resnik_mod_314_log2_v2_2018-11-12.csv",stringsAsFactors = F)) 
n3_100k_max <- read_csv(file.path(base_dir,'100k_3pats_resnik_mod_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n4_100k_max <- read_csv(file.path(base_dir,'100k_4pats_resnik_mod_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n5_100k_max <- read_csv(file.path(base_dir,'100k_5pats_resnik_mod_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))
n7_100k_max <- read_csv(file.path(base_dir,'100k_7pats_resnik_mod_314_log2_v2_2018-11-12.csv',stringsAsFactors = F))

exp314 <- read_csv(file.path(base_dir,'314_Expanded_v2.csv'))
local_IC <- read_csv(file.path(base_dir,"Local_IC_314_log2_v2_2018-11-09.csv")) 

pat_table_base <- read_csv(file.path(base_dir,"Local_Base_per_patient_315_v2_2018-11-09.csv"))
pat_table_prop <- read_csv(file.path(base_dir,"Local_Prop_per_patient_315_v2_2018-11-09.csv"))


path <- read_csv(file.path(base_dir,"paths_data_frame.csv"))
allHPOs <- read_csv(file.path(base_dir,"HPO_is.a_tree.csv"))
hpo_ancs <- read_csv(file.path(base_dir,"HPO_Ancestors_Compressed.csv"))
