######
#Load required libraries
######
library(tidyverse)


######
#STEP 1: Load required files
######
#Freeze files from the 320 trio cohort
setwd("/Volumes/helbig_lab/Dropbox/HPO_collaborated/Permutation_Analysis_Manuscript/expanded_freeze_315_v2")
variant = read.csv("/Volumes/helbig_lab/Dropbox/Freeze/V7_Cumulative_328_Freeze_csv.csv", 
                   stringsAsFactors = F) %>%  dplyr::rename(famID = Family_ID) %>%
  #harmonize EG0192, EG0180 so conistent with other dataframes
  mutate(famID=replace(famID, famID == "EG0192P", "EG0192")) %>% 
  mutate(famID=replace(famID, famID == "EG0180P", "EG0180"))
  
sim_score = read.csv("sim_matrix_resnik_avg_315_log2_2018-11-11.csv",stringsAsFactors = F,row.names = 1)


######
#STEP 2: Clean and organize the data
##Generate required tables and dataframes
######
#Match famIDs in sim_score and variant

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

#Data check
## doublecheck = 320
print(paste("There 314 trios with sim scores?:", 
           ((variant_sim$famID %>% unique %>% length)==314), sep = ""))


######
#STEP 3: Isolate variants
######
# Table of denovo variants
##Setting additional constraints on denovo mutations
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

########
#Function: gene_df
########
#Input: gene of interest from the data set
#Output: return table of famIDs with gene compared to one another (i.e. with their sim scores) 

gene_df <- function(gene)
{
  az <- tab1 %>% filter(Gene.refGene==gene) %>% unique
  len1 = nrow(az)
  
  #matrix of all combinations
  matrix = t(combn(1:len1, 2))
  len2 = nrow(matrix)
  
  df <- as.data.frame(matrix(ncol=4,nrow=len2))
  names(df) <- c('fam1','fam2','gene','sim_score')
  df[,3] <- gene
  
  #fill df with indices of matrix and add sim_score
  for (i in 1:len2)
  {
    x_name = az$famID[matrix[i,1]]
    y_name = az$famID[matrix[i,2]]
    df[i,1] <- x_name
    df[i,2] <- y_name
    y_name_cor <- gsub("-",".",y_name) #column names in dataframe replace "_" with "."
    df[i,4] <- sim_score[x_name,y_name_cor]
  }
  df <- df[!duplicated(df),] #matrix of unique combinations
  return(df)
}

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


#########
#STEP 4: Create similarity score distributions based on patients with the same gene
## Creating large dataframe of pvalues based on similarity between patients with
###the same denovo gene
#########


n1_100k = read.csv("100k_2pats_resnik_avg_315_log2_v2_2018-11-12.csv",stringsAsFactors = F) 

n3_100k = read.csv('100k_3pats_resnik_avg_315_log2_v2_2018-11-12.csv',stringsAsFactors = F)

n4_100k = read.csv('100k_4pats_resnik_avg_315_log2_v2_2018-11-12.csv',stringsAsFactors = F)

n5_100k = read.csv('100k_5pats_resnik_avg_315_log2_v2_2018-11-12.csv',stringsAsFactors = F)

n7_100k = read.csv('100k_7pats_resnik_avg_315_log2_v2_2018-11-12.csv',stringsAsFactors = F)

#Function: sim_random
##Output: picks sim_score from two random pairs
sim_random <- function() {
  random = sample(1:nrow(sim_score),2,replace = FALSE)
  sims_rnd <- sim_score[random[1],random[2]]
  return(sims_rnd)
}

#Function: draw_av
##Output: returns the average sim_score in n random individuals 
draw_av <- function(n) {
  count = 0
  for (i in 1:n) {
    count = count + sim_random()
  }
  count = count/n
  return(count)
}

#Function: draw_median
##Output: returns the median sim_score in n random individuals 
draw_median <- function(n) {
  count = c()
  for (i in 1:n) {
    count[i] = sim_random()
  }
  median_i <- median(count)
  return(median_i)
}

#Function: draw_median
##Output: returns the mode sim_score in n random individuals 
estimate_mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}

#Function: draw_mode
##Output: returns the mode sim_score in n random individuals 
draw_mode <- function(n) {
  count = c()
  for (i in 1:n) {
    count[i] = sim_random()
  }
  mode_i <- estimate_mode(count)
  return(mode_i)
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
#STEP 4: Find the p_value for each gene's similarity
## using the average, median, and mode similarity,for each gene from gene_count table
#########

#Function: exact_p_av2
##Determines p-value for the average
exact_p_av2 <- function(sample_size,sim_score_med,number_random_draws)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  # if (n > 1){
  #   for (i in 1:nx) {
  #     ax[i,1] <- draw_av(n)
  #   }
   #}
  if (n == 3){
    ax <- n3_100k
  }
  if (n == 4){
    ax <- n4_100k
  }
  if (n == 5){
    ax <- n5_100k
  }
  if (n == 7){
    ax <- n7_100k
  }
  if (n == 2){
    ax <- n1_100k
  }
  
  # hist(ax[,1],breaks=40,main= paste("draws n=",n))
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x){
    foo = ecdf(ax[,1])
    return(1 - foo(x))
  }
  
  p <- p_value(sim_score_med)
  return(p)
}

#Function: exact_p_med2
##Determines p-value for the median
exact_p_med2 <- function(sample_size,sim_score_med,number_random_draws)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  if (n == 3){
    ax <- n3_100k
  }
  if (n == 4){
    ax <- n4_100k
  }
  if (n == 5){
    ax <- n5_100k
  }
  if (n == 7){
    ax <- n7_100k
  }
  if (n == 2){
    ax <- n1_100k
  }
  
  # hist(ax[,1],breaks=40,main= paste("draws n=",n))
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {
    foo = ecdf(ax[,1])
    return(1 - foo(x))
  }
  p <- p_value(sim_score_med)
  return(p)
  
}

#Function: exact_p_mod2
##Determines p-value for the mode
exact_p_mod2 <- function(sample_size,sim_score_mod,number_random_draws){
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  nx = number_random_draws
  n = sample_size
  ax <- as.data.frame(matrix(nrow=nx,ncol=1))
  names(ax) <- c("random_draw")
  
  if (n == 3){
    ax <- n3_100k
  }
  if (n == 4){
    ax <- n4_100k
  }
  if (n == 5){
    ax <- n5_100k
  }
  if (n == 7){
    ax <- n7_100k
  }
  if (n == 2){
    ax <- n1_100k
  }
  
  # hist(ax[,1],breaks=40,main= paste("draws n=",n))
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {
    foo = ecdf(ax[,1])
    return(1 - foo(x))
  }
  p <- p_value(sim_score_mod)
  return(p)
}


#loop through each gene for p_av, p_med
for (i in 1:nrow(gene_count)){
  #print(paste(gene_count$gene[i]," - number of pairs",gene_count$pairs[i]))
  e1 = gene_count$n_pats[i] #sample size
  #e2 = gene_count$mode_sim[i] #sim score mod
  e3 = 10000 # 10K unless n=1, then 100K (predefined)
  p_average = exact_p_av2(e1,gene_count$av_sim[i],e3)
  p_median = exact_p_med2(e1,gene_count$median_sim[i],e3)
  p_mod = exact_p_mod2(e1,gene_count$mode_sim[i],e3)
  gene_count[i,7] <- p_average
  gene_count[i,8] <- p_median
  gene_count[i,9] <- p_mod
  print(paste("...p-value average = ",p_average))
  print(paste("...p-value median = ",p_median))
  print(paste("...p-value mode = ",p_mod))
}

#write.table(gene_count,"gene_count_DNM_filtered.txt")
#write.csv(gene_count,"gene_count_res_avg_log2_314_v2_cor.csv",row.names = F)
