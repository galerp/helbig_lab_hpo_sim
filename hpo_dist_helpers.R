
library(tidyverse)
library(memoise)
library(dplyr)


#########
#Function - mica (memoised via memo_mica())
##Find the Most Informative Common Ancestor of 2 HPO terms using Resnik (1995)
##Input - two HPO terms
#         hpo1
#         hpo2
##Output - 
#         mica_ic - the information content of the most informative common ancestor (MICA)
#########

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

pat_base <- function(exp321){
  pat_table_base <- exp321 %>% 
  dplyr::select(famID,HPO) %>% 
separate_rows(HPO, sep = ";") %>% unique()
  return(pat_table_base)
}

pat_prop <- function(pat_table_base){
pat_table_prop <- pat_table_base %>% 
left_join(hpo_ancs %>% dplyr::select(-definition)) %>% 
dplyr::select(famID, Ancestors) %>% 
dplyr::rename(HPO = Ancestors) %>% 
separate_rows(HPO, sep = ";") %>% 
#Remove duplicated HPO terms in each patient
unique
}



#########
#Function - pat_compare
##Find the similarity score between two patients via the sim_max or sim_av (PMID: 16776819) method
##Input - two patient IDs
#         pat1
#         pat2
##Output - 
#         max_complete - the similarity score between the two patients
#########

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


base_IC <- function(pat_table_base) { 
 base_IC <- pat_table_base %>% 
  dplyr::count(HPO) %>% 
  mutate(local.Base.freq = n/length(unique(pat_table_base$famID))) %>% 
  mutate(Base.local.IC = -log2(local.Base.freq)) %>% 
  dplyr::select(-n)
 return(base_IC)
}

prop_IC <- function(pat_table_prop) {
prop_IC <- pat_table_prop %>% 
  dplyr::count(HPO) %>% 
  dplyr::mutate(local.Prop.freq = n/length(unique(pat_table_base$famID))) %>% 
  dplyr::mutate(Propagated.local.IC = -log10(local.Prop.freq)) %>% 
  dplyr::select(-n)
return(prop_IC)
}

local_IC <- function(allHPOs) {
local_IC <- allHPOs %>% 
  dplyr::select(term) %>% 
  full_join(base_IC, by = c('term' = 'HPO')) %>% 
  full_join(prop_IC, by = c('term' = 'HPO')) %>% 
  replace(., is.na(.),0)
return(local_IC)
}
#########
#Function - Compare_Cohort
##Computes the similarity scores between every patient in the cohort
##Input - 
#         cohort_file - every patient in the cohort with their HPO terms in base format
#         
##Output - 
#         sim_score - a n x n similarity matrix with the sim scores of every patient pair in the cohort
#########

Compare_Cohort <- function(cohort_file){
  #create output table
  patients <- unique(cohort_file$famID)
  dimension <- length(unique(patients))
  
  pat_matrix <- as.data.frame(matrix(ncol=dimension, nrow=dimension))
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  pat_matrix <- sapply(1:dimension, function(y) 
    (sapply(1:dimension, function(x) pat_matrix[y,x] <- pat_compare(names(pat_matrix)[y],rownames(pat_matrix)[x]))))
  pat_matrix <- as.data.frame(pat_matrix)
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  return(pat_matrix)
}


#########
#Function - sim_pat_draw
##Randomly select num_pats number of patients N times from sim_score and using their similiarity 
##scores computes the median, mean and mode score among those patients, creating a distribution
##Input - 
#         sim_score - similarity score matrix of all patients in cohort
#         num_pats - number of patients to randomly draw N times
#         N - the number of iterations to randomly draw patients to compos the distributions
##Output - 
#         r_100k - N iterations of median, mean and mode scores among num_pats number of patients
#########

sim_pat_draw <- function(sim_score, num_pats,num_iterations)  {
  
  r_100k = as.data.frame(matrix(nrow = num_iterations, ncol = 3))
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


#########
#Function - gene_df
##Returns a table of the identifiers (famID) and details of all patients with inputted gene
##Input - 
#         gene - gene of interest within the cohort
#         
##Output - 
#         df - a table with all patients with inputted gene and the similarity score between all 
#              said patients
#########

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

denovo <- function(variant_sim){
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
  return(denovo)
}



gene_compute <- function(gene_count){
  
  for (i in 1:nrow(gene_count)) {
    name_x <- subset(pair_corrected,pair_corrected$gene == as.character(gene_count[i,c("gene")]))
    gene_count[i,c('n_pats')] <- combine(name_x$fam1, name_x$fam2) %>% unique %>% length()
    gene_count[i,c('pairs')] <- nrow(name_x)  
    gene_count[i,c('av_sim')] <- sum(name_x$sim_score)/nrow(name_x)
    gene_count[i,c('median_sim')] <- median(name_x$sim_score)
    if (length(name_x$sim_score) > 1) {gene_count[i,c('mode_sim')] <- estimate_mode(name_x$'sim_score')}
    if (length(name_x$sim_score) == 1) {gene_count[i,c('mode_sim')] <- name_x$sim_score}
  }
  
  for (i in 1:nrow(gene_count)) {
    #print(paste(gene_count$gene[i]," - number of pairs",gene_count$pairs[i]))
    e1 <- gene_count$n_pats[i] #sample size
    #e2 = gene_count$mode_sim[i] #sim score mod
    p_average = exact_p(e1, gene_count$av_sim[i])
    p_median = exact_p(e1,gene_count$median_sim[i])
    p_mod = exact_p(e1,gene_count$mode_sim[i])
    gene_count[i,7] <- p_average
    gene_count[i,8] <- p_median
    gene_count[i,9] <- p_mod
    
  }
  return(gene_count)
}
#########
#Function - sim_random
##Randomly selects a patient pair's sim score

##Output - 
#         sims_rnd - a random sim score from sim_score
#########

sim_random <- function() {
  random = sample(1:nrow(sim_score),2,replace = FALSE)
  sims_rnd <- sim_score[random[1],random[2]]
  return(sims_rnd)
}


#########
#Function - estimate_mode
##Estimates the mode of a given vector of number
##Input - 
#         x - a vector of numbers
##Output - 
#         the estimated mode
#########

estimate_mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}


#########
#Function - exact_p_av2
## Calculates the p-value of the gene according the similarity scores of the patients
## with the gene in question and the previously created distribution
##Input - 
#         sample_size - number of patients with the gene in question
#         similarity_score - the 
##Output - 
#         p - p-value based on the mean values
#########

exact_p <- function(sample_size,similarity_score,method)
{
ax <- data.frame()
  if (sample_size == 3){
    ax <- as.data.frame(n3_100k)
  } else if (sample_size == 4){
    ax <- as.data.frame(n4_100k)
  } else  if (sample_size == 5){
    ax <- as.data.frame(n5_100k)
  } else if (sample_size == 7){
    ax <- as.data.frame(n7_100k)
  } else if (sample_size == 2){
    ax <- as.data.frame(n2_100k)
  }
  
  #p_value inner-function to determine the significance of a value of x 
p_value <- function(x, method) {
  
  if (method == "median") {foo <- ecdf(ax[,1])
  } else if (method == "mean"){foo <- ecdf(ax[,2])
  } else if (method == "mode"){foo <- ecdf(ax[,3])
  }
  return(1 - foo(x))
}
  p <- p_value(similarity_score,method)
  return(p)
}
