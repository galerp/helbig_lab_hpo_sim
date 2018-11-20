library(tidyverse)
library(memoise)
library(dplyr)

######
#STEP 1: Create base and prop table
######

pat_table_base <- exp321 %>% 
  dplyr::select(famID,HPO) %>% 
  separate_rows(HPO, sep = ";") %>% unique()

pat_table_prop <- pat_table_base %>% 
  left_join(hpo_ancs %>% dplyr::select(-definition)) %>% 
  dplyr::select(famID, Ancestors) %>% 
  dplyr::rename(HPO = Ancestors) %>% 
  separate_rows(HPO, sep = ";") %>% 
  #Remove duplicated HPO terms in each patient
  unique

######
#STEP 2: Calculate Local IC
######

base_IC <- pat_table_base %>% 
  dplyr::count(HPO) %>% 
  mutate(local.Base.freq = n/length(unique(pat_table_base$famID))) %>% 
  mutate(Base.local.IC = -log2(local.Base.freq)) %>% 
  dplyr::select(-n)

prop_IC <- pat_table_prop %>% 
  dplyr::count(HPO) %>% 
  dplyr::mutate(local.Prop.freq = n/length(unique(pat_table_base$famID))) %>% 
  dplyr::mutate(Propagated.local.IC = -log10(local.Prop.freq)) %>% 
  dplyr::select(-n)

local_IC <- allHPOs %>% 
  dplyr::select(term) %>% 
  full_join(base_IC, by = c('term' = 'HPO')) %>% 
  full_join(prop_IC, by = c('term' = 'HPO')) %>% 
  replace(., is.na(.),0)


write.csv(local_IC, "Local_IC.csv",row.names = F)

#############
#STEP 3: Simlilarity Analysis - sim_max or sim_av
##Input: local_IC
#        pat_table_base
############
ic <- local_IC[,c("term","Propagated.local.IC")]
hpo_all <- allHPOs[,c("term","definition")]
ic2 <- merge(ic,hpo_all,by='term',all.x=TRUE,all.y=FALSE)

names(local_IC)[1] = "HPO"
ic2 <- local_IC[,c("HPO","Propagated.local.IC")]

mica <- function(hpo1, hpo2)
{
  
  
  # path1_unique <- path %>% dplyr::filter(Term == hpo1)
  path1_unique <- path[which(path$Term == hpo1),]
  
  # path2_unique <- path %>% dplyr::filter(Term == hpo2)
  path2_unique <- path[which(path$Term == hpo2),]
  
  joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1') 
  
  joint2 <- joint1 %>% left_join(ic2, by = c('Counter1' = 'HPO')) 
  
  mica_ic <- joint2$Propagated.local.IC %>% max  
  
  
  return(mica_ic)
  
}

memo_mica <- memoise(mica)


pat_compare <- function(pat1, pat2)
{
  # hpo_pat1 <- pat_table_base %>% dplyr::filter(famID == pat1)
  # hpo_pat2 <- pat_table_base %>% dplyr::filter(famID == pat2)
  
  hpo_pat1 <- pat_table_base[which(pat_table_base$famID == pat1),]
  hpo_pat2 <- pat_table_base[which(pat_table_base$famID == pat2),]
  
  #create data frame with HPO of pat1 in x, HPO of pat2 in y
  
  x_length <- length(hpo_pat1$HPO)
  
  y_length <- length(hpo_pat2$HPO)
  
  
  ic_matrix <- as.data.frame(matrix(ncol=x_length, nrow=y_length))
  
  names(ic_matrix) <- hpo_pat1$HPO
  
  rownames(ic_matrix) <- hpo_pat2$HPO
  
  for (i in 1:y_length){
    
    for(j in 1:x_length)
    {
      ic_matrix[i,j] <- memo_mica(hpo_pat2$HPO[i], hpo_pat1$HPO[j])
      
    }
  }
if(input.yaml$algorithm == 1 ){ 
  
  max_col <- apply(ic_matrix,2,max)
  
  max_row <- apply(ic_matrix,1,max)
  
  max_complete <- sum(max_col,max_row)/2
} else{
  
  max_col <- apply(ic_matrix,2,max)
  max_col <- max_col/length(ic_matrix) #KEY DIFFERENCE IN THIS ALGORITHM
  
  max_row <- apply(ic_matrix,1,max)
  max_row <- max_row/nrow(ic_matrix)
  
  
  max_complete <- sum(max_col,max_row)/2
}
  
  
  return(max_complete)
}

############
#Run Similarity Analysis on entire cohort 
############

Compare_Cohort=function(cohort_file){
  #create output table
  patients=unique(cohort_file$famID)
  dimension <- length(unique(patients))
  
  pat_matrix <- as.data.frame(matrix(ncol=dimension, nrow=dimension))
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  pat_matrix = sapply(1:dimension, function(y) 
    (sapply(1:dimension, function(x) pat_matrix[y,x] <- pat_compare(names(pat_matrix)[y],rownames(pat_matrix)[x]))))
  pat_matrix = as.data.frame(pat_matrix)
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  return(pat_matrix)
}


sim_score <- Compare_Cohort(pat_table_base)


write.csv(pat_matrix,file = "sim_matrix.csv",row.names = T)


###########
#STEP 4: Create Distribution of Similarity Scores
##Input:
#       sim_score
##Output:
#       N_100k - 100k distributions of N patients 
###########


rownames(sim_score) = names(sim_score)


estimate_mode <- function(x){  
  d <- density(x)
  d$x[which.max(d$y)]
}


sim_pat_draw = function(sim_score, num_pats)  {
  
  r_100k = as.data.frame(matrix(nrow = 100000, ncol = 3))
  names(r_100k) = c("median","mean", "mode")
  
  pat_vect = names(sim_score)
  
  for(n in 1: nrow(r_100k)){
    IDs = sample(pat_vect, num_pats)
    sub_sim = sim_score[(rownames(sim_score) %in% IDs), (names(sim_score) %in% IDs)]
    diag(sub_sim) = 12345
    vect_scores =  unlist(sub_sim)
    vect_scores = vect_scores[-which(vect_scores == 12345)]
    
    r_100k$median[n] = median(vect_scores)
    r_100k$mean[n] = mean(vect_scores)
    r_100k$mode[n] = estimate_mode(vect_scores)
  }
  return(r_100k)
}

pat_pair = c(2,3,4,5,7)
for (i in 1:length(pat_pair)) { 
  paste0("n",pat_pair[i],"_100k") = sim_pat_draw(sim_score, pat_pair[i])
 write.csv(paste0("n",pat_pair[i],"_100k"), paste0("100k_",pat_pair[i],"pats_resnik_mod.csv", row.names = F))
}



###########
#STEP 4: Generate P-values for Genes
##Input:
#       N_100k - distributions
#       variant - de novo variants in cohort
##Output:
#       gene_count
###########

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
#Find the p_value for each gene's similarity
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
    ax <- n2_100k
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
    ax <- n2_100k
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
    ax <- n2_100k
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

}


write.csv(gene_count,"gene_count.csv",row.names = F)
