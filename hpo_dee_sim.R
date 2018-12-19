
library(tidyverse)
library(memoise)

if(file.exists("input.yml") == T){
  input.yaml <- yaml::read_yaml("input.yml")
}else{
  message('Input YAML not found. Please Save the TAML file in current Directory\n')
  break;
}

if(is.null(input.yaml$patient_phenome) == F){
  exp321 <- read_csv(input.yaml$patient_phenome)
}else{
  message(' Please mention the Patient Phenome File in the specified format - Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$hpo_tree) == F){
  allHPOs <- read_csv(input.yaml$hpo_tree)
}else{
  message(' Please mention the HPO tree File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$hpo_ancestor) == F){
  hpo_ancs <- read_csv(input.yaml$hpo_ancestor)
}else{
  message(' Please mention the HPO Ancestor File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}

options(stringsAsFactors = F)


if(is.null(input.yaml$hpo_path) == F){
  path <- read_csv(input.yaml$hpo_path)
}else{
  message(' Please mention the HPO Path File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$variant_file) == F){
  variant <- read_csv(input.yaml$variant_file)
}else{
  message(' Please mention the Variants File in the specified format - Cant Proceed without that \n')
  break;
}



if( is.null(input.yaml$local_IC)) { 
  source("hpo_dist.R")
} else if(is.null(input.yaml$sim_matrix)) {
  
  local_IC <- read_csv(input.yaml$local_IC)
  pat_table_base <- read_csv(input.yaml$patient_hpo_base)
  pat_table_prop <- read_csv(input.yaml$patient_hpo_prop)
  source("precalc_IC.R")
} else {
  local_IC <- read_csv(input.yaml$local_IC)
  pat_table_base <- read_csv(input.yaml$patient_hpo_base)
  pat_table_prop <- read_csv(input.yaml$patient_hpo_prop)
  sim_score <- read.csv(input.yaml$sim_matrix, row.names = 1)
  source("precalc_sim_matrix.R")
}


