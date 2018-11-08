#Resnik default



#Different MICA algorithm experimentation
library(memoise)
library(dplyr)



#Add Data
exp321=read.csv("320_Expanded_v1_PG_clean.csv",stringsAsFactors = F)

#allHPOs - static file
allHPOs=read.csv("HPO_is.a_tree.csv",stringsAsFactors = F)
allHPOs=allHPOs[!duplicated(allHPOs$term),]

#local_IC - use most recent version 
local_IC=read.csv("Local_IC_320_v1_log10.csv",stringsAsFactors = F) 


#HPO terms in cohort (base and propagated)
pat_table_base <- read.csv("Local_Base_per_patient_320_v2.csv")
pat_table_prop <- read.csv("Local_Prop_per_patient_320_v1.csv")

#generate HPO/IC table (required for mica function)
ic <- local_IC[,c("term","Propagated.local.IC")]
hpo_all <- allHPOs[,c("term","definition")]
ic2 <- merge(ic,hpo_all,by='term',all.x=TRUE,all.y=FALSE)


path = read.csv("paths_data_frame.csv",stringsAsFactors = F)


ancest_comp=read.csv("HPO_Ancestors_Compressed.csv",stringsAsFactors = F)

#generate HPO/IC table (required for mica function)
names(local_IC)[1] = "HPO"
ic2 <- local_IC[,c("HPO","Propagated.local.IC")]



################ 
# MICA function
################ 
mica <- function(hpo1, hpo2)
{
  
  
 # path1_unique <- path %>% dplyr::filter(Term == hpo1)
   path1_unique <- path[which(path$Term == hpo1),]
  
 # path2_unique <- path %>% dplyr::filter(Term == hpo2)
   path2_unique <- path[which(path$Term == hpo2),]
  
  joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1') #merge(path1_unique,path2_unique,by="Counter1",all=F)
  
  joint2 <- joint1 %>% left_join(ic2, by = c('Counter1' = 'HPO')) #merge(joint1,ic2,by.x='Counter1',by.y='HPO',all.x=TRUE,all.y=FALSE)

  mica_ic <- joint2$Propagated.local.IC %>% max  #max(joint2$prop.IC)
  
  
  return(mica_ic)
  
}

memo_mica <- memoise(mica)



################ 
#Patient Comparison Function
################ 

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
      # HP = mixedsort(hpo_pat2$HPO[i],hpo_pat1$HPO[j])  #this will make the comparison memoise faster
      ic_matrix[i,j] <- memo_mica(hpo_pat2$HPO[i], hpo_pat1$HPO[j])
      
    }
  }
  
  
  max_col <- apply(ic_matrix,2,max)
  
  max_row <- apply(ic_matrix,1,max)
  
  max_complete <- sum(max_col,max_row)/2
  
  return(max_complete)
}


############### run on entire cohort 
Compare_Cohort=function(cohort_file){
  #create output table
  patients=unique(cohort_file$famID)
  dimension <- length(unique(patients))

  pat_matrix <- as.data.frame(matrix(ncol=dimension, nrow=dimension))
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  pat_matrix = sapply(1:dimension, function(y) (sapply(1:dimension, function(x) pat_matrix[y,x] <- pat_compare(names(pat_matrix)[y],rownames(pat_matrix)[x]))))
  pat_matrix = as.data.frame(pat_matrix)
  
  names(pat_matrix) <- patients[1:dimension]
  rownames(pat_matrix) <- patients[1:dimension]
  
  return(pat_matrix)
}


pat_matrix=Compare_Cohort(pat_table_base)

write.csv(pat_matrix,file = "full_trio_sim_matrix_resnik_def_fast_2018-10-30.csv",row.names = T)
