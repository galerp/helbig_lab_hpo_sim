
library(tidyverse)
library(memoise)

capture <- commandArgs(trailingOnly = TRUE)

opt1 = list(make_option(c("--input"), type = "character", default = "input.yml", dest = "input"))


user_input <- function(name, argv) {
  return(name %in% names(argv))
}

argv <- parse_args(OptionParser(option_list = opt1))

if (user_input("input", argv)) {
  input = argv$input 
    if(file.exists(input) == T){
      input.yaml <- yaml::read_yaml(input)
    }else{
      message('\n Input YAML not found \n')
      break;
    }
  } else {
  message('Cannot proceed without input yaml file. Please use "--input" flag .\n')
}


if(is.null(input.yaml$output_dir) == T){
  message('\n Please mention the Field output_dir in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$patient_phenome) == F){
  exp321 <- read_csv(input.yaml$patient_phenome)
}else{
  message('\n Please mention the Patient Phenome File in the specified format - Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$hpo_tree) == F){
  allHPOs <- read_csv(input.yaml$hpo_tree)
}else{
  message('\n  Please mention the HPO tree File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$hpo_ancestor) == F){
  hpo_ancs <- read_csv(input.yaml$hpo_ancestor)
}else{
  message('\n  Please mention the HPO Ancestor File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}

options(stringsAsFactors = F)


if(is.null(input.yaml$hpo_path) == F){
  path <- read_csv(input.yaml$hpo_path)
}else{
  message('\n  Please mention the HPO Path File (Should be able to get it from here "link to be added")- Cant Proceed without that \n')
  break;
}


if(is.null(input.yaml$variant_file) == F){
  variant <- read_csv(input.yaml$variant_file)
}else{
  message('\n  Please mention the Variants File in the specified format - Cant Proceed without that \n')
  break;
}


if( is.null(input.yaml$local_IC)) { 
 message("\n  The minimum required input files were read. The program will run entirely, generating these files : computing the Information Content, Similarity Matrix and Permutation Analysis \n  \n ")
  source("hpo_dist.R")
} else if(is.null(input.yaml$sim_matrix)) {
  
message("\n  The minimum required input files along with the Local_IC, hpo_base and hpo_prop  were read. The program will only run the steps for  generating these files : computing Similarity Matrix and Permutation Analysis \n  \n  ")
  local_IC <- read_csv(input.yaml$local_IC)
  pat_table_base <- read_csv(input.yaml$patient_hpo_base)
  pat_table_prop <- read_csv(input.yaml$patient_hpo_prop)
  source("precalc_IC.R")
} else {
 
message("\n  The minimum required input files along with the Local_IC, hpo_base, hpo_prop and sim_matrix  were read. The program will only run the steps for  generating the files from Permutation Analysis \n  \n ")
  local_IC <- read_csv(input.yaml$local_IC)
  pat_table_base <- read_csv(input.yaml$patient_hpo_base)
  pat_table_prop <- read_csv(input.yaml$patient_hpo_prop)
  sim_score <- read.csv(input.yaml$sim_matrix, row.names = 1)
  source("precalc_sim_matrix.R")
}


