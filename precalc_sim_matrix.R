
source("hpo_dist_helpers.R")
rownames(sim_score) = names(sim_score)

# pat_pair = c(2,3,4,5,7)
# for (i in 1:length(pat_pair)) { 
#   paste0("n",pat_pair[i],"_100k") = sim_pat_draw(sim_score, pat_pair[i])
#   write.csv(paste0("n",pat_pair[i],"_100k"), paste0("100k_",pat_pair[i],"pats_resnik_mod.csv", row.names = F))
# }
num_iterations = input.yaml$num_of_iterations

n2_100k <-  sim_pat_draw(sim_score, 2,num_iterations) # read_csv("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_314_v2/100k_2pats_resnik_avg_314_log2_v2_2018-11-12.csv") 
n3_100k <- sim_pat_draw(sim_score, 3,num_iterations) # read_csv("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_314_v2/100k_3pats_resnik_avg_314_log2_v2_2018-11-12.csv") 
n4_100k <- sim_pat_draw(sim_score, 4,num_iterations) # read_csv("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_314_v2/100k_4pats_resnik_avg_314_log2_v2_2018-11-12.csv") 
n5_100k <- sim_pat_draw(sim_score, 5,num_iterations) # read_csv("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_314_v2/100k_5pats_resnik_avg_314_log2_v2_2018-11-12.csv") 
n7_100k <- sim_pat_draw(sim_score, 7,num_iterations) # read_csv("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_314_v2/100k_7pats_resnik_avg_314_log2_v2_2018-11-12.csv") 


###########
#STEP 4: Generate P-values for Genes
##Input:
#       N_100k - distributions
#       variant - de novo variants in cohort
##Output:
#       gene_count
###########
names(variant)[1] <- "famID"
famIDs_var <- variant$famID %>% unique 
famIDs_var <- famIDs_var %>% as.data.frame %>% 
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


denovo <- variant_sim %>%  
  filter(AD2_Proband >=10 & 
           0.75 >= Ratio_Proband &
           Ratio_Proband >= 0.25 &
           Probability <= 70 &
           Type == "Denovo" & 
           is.na(esp6500siv2_all) & 
           is.na(genomicSuperDups) & 
           is.na(X1000g2015aug_all) & 
           is.na(ExAC_ALL) ) 

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
  print(i)
  print(all_genes$gene[i])
  temp <- gene_df(all_genes$gene[i])
  pair_corrected <- rbind(pair_corrected,temp)
}



#Determine number of pairs per gene
gene_x <- unique(pair_corrected$gene)

gene_count <- as.data.frame(matrix(ncol=9,nrow=length(gene_x))) 

names(gene_count) <- c("gene","n_pats","pairs","av_sim","median_sim","mode_sim","p_av","p_median","p_mode")

gene_count <- gene_count %>% mutate(gene = gene_x)


#Finding the average, median, and mode similarity,for each gene in the cohort
for (i in 1:nrow(gene_count)) {
  name_x <- subset(pair_corrected,pair_corrected$gene == as.character(gene_count[i,c("gene")]))
  gene_count[i,c('n_pats')] <- combine(name_x$fam1, name_x$fam2) %>% unique %>% length()
  gene_count[i,c('pairs')] <- nrow(name_x)  
  gene_count[i,c('av_sim')] <- sum(name_x$sim_score)/nrow(name_x)
  gene_count[i,c('median_sim')] <- median(name_x$sim_score)
  if (length(name_x$sim_score) > 1) {gene_count[i,c('mode_sim')] <- estimate_mode(name_x$'sim_score')}
  if (length(name_x$sim_score) == 1) {gene_count[i,c('mode_sim')] <- name_x$sim_score}
}


#########
#Find the p_value for each gene's similarity
## using the average, median, and mode similarity,for each gene from gene_count table
#########


#loop through each gene for p_av, p_med
for (i in 1:nrow(gene_count)) {
  #print(paste(gene_count$gene[i]," - number of pairs",gene_count$pairs[i]))
  e1 <- gene_count$n_pats[i] #sample size
  #e2 = gene_count$mode_sim[i] #sim score mod
  p_average <-  exact_p(e1, gene_count$av_sim[i])
  p_median <- exact_p(e1,gene_count$median_sim[i])
  p_mod <- exact_p(e1,gene_count$mode_sim[i])
  gene_count[i,7] <- p_average
  gene_count[i,8] <- p_median
  gene_count[i,9] <- p_mod
  
}


write.csv(gene_count,"gene_count_1.csv",row.names = F)
