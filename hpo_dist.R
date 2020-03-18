library(tidyverse)
library(memoise)

start <- Sys.time()
source("hpo_dist_helpers.R")
######
#STEP 1: Create base and prop table
######
message(" \n Step 1 is executed: HPO Base and Prop Table is created \n \n ")

pat_table_base <- pat_base(patient_phenome) 
pat_table_prop <- pat_prop(pat_table_base)

######
#STEP 2: Calculate Local IC
######
message("\n  Step 1 is done. Calculating Local Information Content \n  \n ")

base_IC <- base_calc_IC(pat_table_base)
prop_IC <- prop_calc_IC(pat_table_prop)
local_IC <- local_calc_IC(allHPOs)

write.csv(local_IC, paste0(input.yaml$output_dir,"/Local_IC.csv"),row.names = F)

#############
#STEP 3: Simlilarity Analysis - sim_max or sim_cm
############

message("\n  Local_IC file is written. Calculating the similarity_matrix \n  \n ")

ic <- local_IC[,c("term","Propagated.local.IC")]
hpo_all <- allHPOs[,c("term","definition")]
ic2 <- merge(ic,hpo_all,by='term',all.x=TRUE,all.y=FALSE)

names(local_IC)[1] = "HPO"
ic2 <- local_IC[,c("HPO","Propagated.local.IC")] 


memo_mica <- memoise(mica)

############
#Run Similarity Analysis on entire cohort 
############

sim_score <- Compare_Cohort(pat_table_base)


write.csv(sim_score,paste0(input.yaml$output_dir,"/sim_matrix.csv"),row.names = T)


###########
#STEP 4: Create Distribution of Similarity Scores
###########

message("\n  Sim score file is written. Permutation Analysis is to be followed \n  \n ")
famIDs_var <- variant$famID %>% unique %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>% dplyr::mutate(var = "variant")

famIDs_sim <- sim_score %>% rownames %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>%  dplyr::mutate(sim = "sim_score")

fam_combined <- famIDs_sim %>%  dplyr::full_join(famIDs_var)


#Filtering data
##Trios without sim_scores and without variants
no_sim <- as.vector(fam_combined$famID[is.na(fam_combined$sim) == TRUE])
no_var <- as.vector(fam_combined$famID[is.na(fam_combined$var) == TRUE])

if(length(no_sim) !=0){
  message(paste0("\n  The following IDs are listed in the variant file, but do not have a sim score: \n",
                no_sim, " \n  \n "))
  }

message('IDs from sim matrix who do not have a variant : \n')
  message( paste(rownames(sim_score)[rownames(sim_score) %!in% famIDs_var$famID][1:10],"\n"), '[truncated]...', length(rownames(sim_score)[rownames(sim_score) %!in% famIDs_var$famID]) - 10,'more') #sim_score but no variant



#Trios with sim_scores
all_sim <- as.vector(fam_combined$famID[is.na(fam_combined$sim) == FALSE])
variant_sim <- variants %>% filter(famID %in% all_sim)


#Table of denovos with famID and gene
tab1 <- variant_sim %>% dplyr::select(famID, Gene.refGeneWithVer) %>% unique


#list of all genes
all_genes <- tab1 %>% dplyr::count(Gene.refGene) %>% 
  dplyr::rename(gene = Gene.refGene, Freq = n) %>% 
  filter(Freq > 1)


num_iterations = input.yaml$num_of_iterations

variants = input.yaml$variant_file

# Calculate the number of cores available
n_cores <- detectCores() - 1

sim_pat_perm <- sim_par_perm(sim_score, all_genes, n_cores, N)

saveRDS(sim_pat_perm, paste0(input.yaml$output_dir,"100k_sim_permutation.rds"))

###########
#STEP 5: Generate P-values for Semantic Similarity within
###########

message("\n Permutations are done.Saving output in '100k_sim_permutation.rds' \n
          p-value for gene to phenotype are being calculating \n \n  ")


#Creating dataframe of similarity comparisons with every combination of patient pairs 
##with the same denove gene

#initialize with first gene
pair_corrected <- gene_df(all_genes$gene[1])

#Create for all genes
for (i in 2:nrow(all_genes)){
  temp <- gene_df(all_genes$gene[i])
  pair_corrected <- pair_corrected %>% bind_rows(temp)
}



#Determine number of pairs per gene
gene_x <- unique(pair_corrected$gene)

gene_count <- as.data.frame(matrix(ncol=9,nrow=length(gene_x))) 

names(gene_count) <- c("gene","n_pats","pairs","av_sim","median_sim","mode_sim","p_av","p_median","p_mode")

gene_count <- gene_count %>% mutate(gene = gene_x)


#Finding the average, median, and mode similarity,for each gene in the cohort
#########
#Find the p_value for each gene's similarity
## using the average, median, and mode similarity,for each gene from gene_count table
#########

#loop through each gene for p_av, p_med
gene_stat = gene_compute(gene_count) 


write.csv(gene_stat,paste0(input.yaml$output_dir,"/gene_count.csv"),row.names = F)


message("\n p-value for gene to phenotype has been calculated and saved in file "gene_count.csv \n \n  ")

###########
#STEP 5: Generate P-values for Genes
###########


message("\n Generating p-value for gene using the Denovolyzer \n \n  ")




message("\n  The Entire script is now run successfully. Please find the Final output file gene_count.csv in the output directory \n ")
stop = Sys.time()
stop - start

