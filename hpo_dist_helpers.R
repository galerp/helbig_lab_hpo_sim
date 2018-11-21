
library(tidyverse)
library(memoise)
library(dplyr)

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


sim_random <- function() {
  random = sample(1:nrow(sim_score),2,replace = FALSE)
  sims_rnd <- sim_score[random[1],random[2]]
  return(sims_rnd)
}


#Function: draw_median
##Output: returns the mode sim_score in n random individuals 
estimate_mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}


 #Functionº : exact_p_av2
##Determines p-value for the average
exact_p <- function(sample_size,similarity_score)
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
  
  # hist(ax[,1],breaks=40,main= paste("draws n=",n))
  
  #p_value inner-function to determine the significance of a value of x 
p_value <- function(x) {
  foo = ecdf(ax[,1])
  return(1 - foo(x))
}
  p <- p_value(similarity_score)
  return(p)
}
