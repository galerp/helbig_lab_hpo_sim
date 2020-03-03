######
#Load required libraries
######
# library(clustermq)
library(dplyr)
library(doParallel)
source("sim_pat_draw.R")
#library(Hmisc)
'%!in%' <- function(x,y)!('%in%'(x,y))
#USACHP180495278180


######
#STEP 1: Load required files
######
#Big Chill dataset

# setwd("/mnt/isilon/helbig_lab/Dropbox/galerp/HPO_testing/EGRP_EPGP_sim/EGRP_EPGP_v10")
variants = read.csv("/mnt/isilon/helbig_lab/Users/ganesans/Trio_Analysis/Reannotate/846_filtered_denovo_genes.csv",stringsAsFactors = F) %>%
  rename(famID = Pat) %>% 
  # separate_rows(sep = ";", famID) %>% 
  mutate(famID = gsub("Fires","FIRES",famID)) %>%
  mutate(famID = gsub("FIRES48","FIRES_48",famID)) %>%
  mutate(famID = gsub("FIRES_32","FIRES32",famID)) %>%
  mutate(famID = gsub("KIEL_TRIO","Kiel_Trio",famID)) %>%
  mutate(famID = gsub("ISR_172","ISR172",famID)) %>% 
  mutate(famID = gsub("-",".",famID)) 

sim_score = read.csv("/mnt/isilon/helbig_lab/Dropbox/galerp/HPO_testing/EGRP_EPGP_sim/EGRP_EPGP_v10/mod_sim_epgp_egrp_v10.csv",stringsAsFactors = F) %>% select(-X)
sim_score<- do.call(data.frame, lapply(sim_score, function(x) {
  replace(x, is.infinite(x) | is.na(x), 0)
})
)
diag(sim_score) <- 0
cor_names = names(sim_score)
names(sim_score) = cor_names
rownames(sim_score) = cor_names
rm(cor_names)
###########################################################

######
#STEP 2: Clean and organize the data
##Generate required tables and dataframes
######
#Unpack fill with famIDs and their corresponding variants
famIDs_var <- data.frame(variant = rep(variants$Gene.refGeneWithVer, sapply(strsplit(variants$famID, split = ";"), length)), famID = unlist(strsplit(variants$famID, split = ";")))
famIDs_var <- famIDs_var[!duplicated(famIDs_var),]

famIDs_sim <- sim_score %>% rownames %>% as.data.frame %>% 
  dplyr::rename('famID' = '.') %>%  dplyr::mutate(sim = "sim_score")

fam_combined <- famIDs_sim %>%  dplyr::left_join(famIDs_var) %>% 
  filter(!is.na(sim), !is.na(variant)) %>% #only famIDs in the sim matrix (should be all)
  select(famID, sim, variant) %>% 
  unique #remove instances where patient has same gene which results in duplications



#DATA CHECK
##Trios without sim_scores or without variants

message('Checking if all famID from vcf are present in sim matrix \n')
if( length(famIDs_var$famID[famIDs_var$famID %!in% rownames(sim_score)]) != 0 ){
  message('sim matrix  do not have these famID that are present in vcf\n')
  cat(famIDs_var$famID[famIDs_var$famID %!in% rownames(sim_score)])
  break;
}#variant but no sim_score

message('IDs from sim matrix who do not have a variant : \n')
  message( paste(rownames(sim_score)[rownames(sim_score) %!in% famIDs_var$famID][1:10],"\n"), '[truncated]...', length(rownames(sim_score)[rownames(sim_score) %!in% famIDs_var$famID]) - 10,'more') #sim_score but no variant



######
#STEP 3: Isolate variants
######
# Table of denovo variants

#list of all genes that occur in more than one patient 
all_genes <- fam_combined %>%  dplyr::count(variant) %>% 
  dplyr::rename(Freq = n) #%>% 
# filter(Freq > 1) %>% 
# filter(Freq < 25)


#Create all 100k files- sending each in as a seperate job in the cluster
pat_pairs = all_genes$Freq %>% unique %>% as.list


# for( i in 1:length(pat_pairs)) {
#   n = pat_pairs[i]
#   destfile = paste("100k_",n,"mod_epgp-egrp_v10.csv", sep = "")
#   
#   if(!file.exists(destfile)){
#     print(n)
#     ndraw = paste("n100k = sim_pat_draw(sim_score,", n, "); 
#                   write.csv(n100k,","'",destfile,"',","row.names = F)",sep="")
#     rname = paste("mod_100k_",n,".R",sep = "")
#     p100kR = paste("source('pat_100k_draw_npat_res.R');",ndraw)
#     write(p100kR,rname)
# 
#     shfile = paste0("/cm/shared/apps_chop/R/3.6.0/bin/Rscript ", rname)
#     shname = paste("p",n,"_100k_mod.sh",sep = "")
#     write(shfile,shname)
#     # system(paste0("qsub -cwd -l mem_free=5g,h_vmem=5g -V ",shname))
#   }}

# options(
#       clustermq.template = paste0("sge_template")
#     )
# 
# foo <- "/mnt/isilon/helbig_lab/Dropbox/galerp/HPO_testing/EGRP_EPGP_sim/EGRP_EPGP_v10/mod_sim_epgp_egrp_v10.csv"
# 
# Q(sim_pat_draw,pkgs="dplyr",task_id=1:length(pat_pairs),n_jobs=length(pat_pairs), max_calls_worker = 1, template=list(job_name="compute_coverage"), export = list(pat_pairs = pat_pairs), sim_score_file =foo)

registerDoParallel(makeCluster(4))
start_time <- Sys.time()
## creates the 100k dataframe for each pat_pairs and stores as a list
mod_sim_pat <- foreach( j = 1:length(pat_pairs)) %dopar% { 
  library(dplyr)
  sim_pat_draw(sim_score_file = "/mnt/isilon/helbig_lab/Dropbox/galerp/HPO_testing/EGRP_EPGP_sim/EGRP_EPGP_v10/mod_sim_epgp_egrp_v10.csv",task_id = j, permutation = 100000)
}
end_time <-  Sys.time()
end_time - start_time
names(mod_sim_pat) <- paste0("n_pats_",unlist(pat_pairs))
saveRDS(mod_sim_pat,"100k_mod_epgp-egrp_v10.rds")




########
#Function: gene_df
########
#Input: gene of interest from the data set
#Output: return table of famIDs with gene compared to one another (i.e. with their sim scores) 

gene_df <- function(gene)
{
  az <- fam_combined %>% filter(variant==gene) %>% unique
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
    #y_name_cor <- gsub("-",".",y_name) #column names in dataframe replace "_" with "."
    df[i,4] <- sim_score[x_name,y_name]
    if(is.na(df[i,4])){
      print(x_name)
      print(y_name)
    }
  }
  df <- df[!duplicated(df),] #matrix of unique combinations
  return(df)
}

#Creating dataframe of similarity comparisons with every combination of patient pairs 
##with the same denove gene

#initialize with first gene
pair_corrected <- gene_df(all_genes$variant[1])

#Create for all genes
for (i in 2:nrow(all_genes)){

  temp <- gene_df(all_genes$variant[i])
  pair_corrected <- rbind(pair_corrected,temp)
}
pair_corrected <- pair_corrected %>% filter(!is.na(sim_score))

#########
#STEP 4: Create similarity score distributions based on patients with the same gene
## Creating large dataframe of pvalues based on similarity between patients with
###the same denovo gene
#######




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
  gene_count[i,c('n_pats')] <- dplyr::combine(name_x$fam1, name_x$fam2) %>% unique %>% length()
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
exact_p_av2 <- function(sample_size,sim_score_med)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x){
    foo = ecdf(mod_sim_pat[[paste0("n_pats_",sample_size)]]$mean)
    return(1 - foo(x))
  }
  
  p <- p_value(sim_score_med)
  return(p)
}

#Function: exact_p_med2
##Determines p-value for the median
exact_p_med2 <- function(sample_size,sim_score_med)
{
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {
    foo = ecdf(mod_sim_pat[[paste0("n_pats_",sample_size)]]$median)
    return(1 - foo(x))
  }
  p <- p_value(sim_score_med)
  return(p)
  
}

#Function: exact_p_mod2
##Determines p-value for the mode
exact_p_mod2 <- function(sample_size,sim_score_mod){
  #builds matrix with nx random draws with n random individuals, example SCN1A with 21 pairs
  
  #p_value inner-function to determine the significance of a value of x 
  p_value <- function(x) {
    foo = ecdf(mod_sim_pat[[paste0("n_pats_",sample_size)]]$mode)
    return(1 - foo(x))
  }
  p <- p_value(sim_score_mod)
  return(p)
}


#loop through each gene for p_av, p_med
for (i in 1:nrow(gene_count)){

  e1 = gene_count$n_pats[i] #sample size
  #e2 = gene_count$mode_sim[i] #sim score mod
  # e3 = 100000 # 10K unless n=1, then 100K (predefined)
  p_average = exact_p_av2(e1,gene_count$av_sim[i])
  p_median = exact_p_med2(e1,gene_count$median_sim[i])
  p_mod = exact_p_mod2(e1,gene_count$mode_sim[i])
  gene_count[i,7] <- p_average
  gene_count[i,8] <- p_median
  gene_count[i,9] <- p_mod

}

 write.csv(gene_count,"gene_count_mod_epgp_egrp_v10.csv",row.names = F)

