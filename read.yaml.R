
library(tidyverse)
library(memoise)

if(file.exists("input.yml") == T){
  input.yaml <- yaml::read_yaml("input.yml")
}else{
  message('Input YAML not found.\n')
  break;
}
options(stringsAsFactors = F)

setwd(input.yaml$working_dir)
exp321 <- read_csv(input.yaml$patient_phenome)
allHPOs <- read_csv(input.yaml$hpo_tree)
hpo_ancs <- read_csv(input.yaml$hpo_ancestor)

  
path <- read_csv(input.yaml$hpo_path)

variant <- read_csv(input.yaml$variant_file)


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


