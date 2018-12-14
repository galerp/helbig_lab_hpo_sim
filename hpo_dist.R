library(tidyverse)
library(memoise)

start <- Sys.time()
source("hpo_dist_helpers.R")
######
#STEP 1: Create base and prop table
######

exp321 <- exp321 %>% mutate(famID = gsub("-","_", famID))

pat_table_base <- pat_base(exp321) 
pat_table_prop <- pat_prop(pat_table_base)

######
#STEP 2: Calculate Local IC
######

base_IC <- base_calc_IC(pat_table_base)
prop_IC <- prop_calc_IC(pat_table_prop)
local_IC <- local_calc_IC(allHPOs)

write.csv(local_IC, "Local_IC.csv",row.names = F)

#############
#STEP 3: Simlilarity Analysis - sim_max or sim_av
############
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


write.csv(sim_score,"sim_matrix.csv",row.names = T)


###########
#STEP 4: Create Distribution of Similarity Scores

###########

rownames(sim_score) = names(sim_score)


num_iterations = input.yaml$num_of_iterations

n2_100k <- sim_pat_draw(sim_score, 2,num_iterations)
n3_100k <- sim_pat_draw(sim_score, 3,num_iterations)
n4_100k <- sim_pat_draw(sim_score, 4,num_iterations)
n5_100k <- sim_pat_draw(sim_score, 5,num_iterations)
n7_100k <- sim_pat_draw(sim_score, 7,num_iterations)


###########
#STEP 4: Generate P-values for Genes

###########

names(variant)[1] <- "famID"
variant <- variant %>% mutate(famID = gsub("-","_", famID))

famIDs_var <- variant$famID %>% unique  %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>% dplyr::mutate(var = "variant")

famIDs_sim <- sim_score %>% rownames %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>%  dplyr::mutate(sim = "sim_score")

fam_combined <- famIDs_sim %>%  dplyr::full_join(famIDs_var)


#Filtering data
##Trios without sim_scores and without variants
no_sim <- as.vector(fam_combined$famID[is.na(fam_combined$sim) == TRUE])
no_var <- as.vector(fam_combined$famID[is.na(fam_combined$var) == TRUE])
#Trios with sim_scores
all_sim <- as.vector(fam_combined$famID[is.na(fam_combined$sim) == FALSE])
variant_sim <- variant %>% filter(famID %in% all_sim)


denovo <- denovo_calc(variant_sim)

#Table of denovos with famID and gene
tab1 <- denovo %>%  dplyr::select(famID, Gene.refGene) %>% unique


#list of all genes
all_genes <- tab1 %>%  dplyr::count(Gene.refGene) %>% 
  dplyr::rename(gene = Gene.refGene, Freq = n) %>% 
  filter(Freq > 1)


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


write.csv(gene_stat,"gene_count.csv",row.names = F)

stop = Sys.time()
stop - start

