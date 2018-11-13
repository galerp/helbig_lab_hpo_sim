# assess pat_compare function
# two different functions
# resnik_max and resnik_average

# resnik_max = "default", "mod" in prior scripts

#use two patients with AP2M1, MAE4, MAE63

#use script from pat_compare function

pat1 = "MAE4"
pat2 = "MAE63"

hpo_pat1 <- pat_table_base[which(pat_table_base$famID == pat1),]
hpo_pat2 <- pat_table_base[which(pat_table_base$famID == pat2),]


#what are the terms for pat1?

hpo_pat1 %>% left_join(allHPOs, by=c('HPO'='term')) -> pat1_def 
pat1_def %>% left_join(local_IC, by='HPO') -> pat1_def2
rm(pat1_def)

#are there duplicate terms for pat1 ?
duplicated(pat1_def2$HPO)

#but parent terms can be duplicated for pat1 
duplicated(pat1_def2$is.a)

# !!! test routine - what are the terms for pat2?

hpo_pat2 %>% left_join(allHPOs, by=c('HPO'='term')) -> pat2_def 
pat2_def %>% left_join(local_IC, by='HPO') -> pat2_def2
rm(pat2_def)

#are there duplicate terms for pat1 ?
duplicated(pat2_def2$HPO)

#but parent terms can be duplicated for pat1 
duplicated(pat2_def2$is.a)

# back to pat_compare

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
#ic_matrix is generate 
# !!! test routine 

#does number of columns, rows fit with the HPO terms per patient

#number of HPO terms in pat1 = "MAE4"
NROW(hpo_pat1)

#number of HPO terms in pat2 = "MAE63"
NROW(hpo_pat2)

#ic_matrix
dim(ic_matrix)[1] == NROW(hpo_pat2) #pat2 in rows
dim(ic_matrix)[2] == NROW(hpo_pat1) #pat1 in columns

#check three individual HPO combinations, example 1
# HP:0001290 in pat2 (rows) and HP:0002342 in pat1 (columns)

#Are both terms present in pat1 and pat2?
"HP:0001290" %in% hpo_pat2$HPO
"HP:0002342" %in% hpo_pat1$HPO

mica("HP:0001290","HP:0002342")
ic_matrix["HP:0001290","HP:0002342"]
subset(allHPOs,allHPOs$term == "HP:0001290")$definition
subset(allHPOs,allHPOs$term == "HP:0002342")$definition
#both meet at phenotypic abnormality

#check three individual HPO combinations, example 2

#Are both terms present in pat1 and pat2?
"HP:0000729" %in% hpo_pat2$HPO
"HP:0000718" %in% hpo_pat1$HPO

mica("HP:0000729","HP:0000718")
ic_matrix["HP:0000729","HP:0000718"]
subset(allHPOs,allHPOs$term == "HP:0000729")$definition
subset(allHPOs,allHPOs$term == "HP:0000718")$definition

#both meet at behavioral abnormality 

#check three individual HPO combinations, example 3

#Are both terms present in pat1 and pat2?
"HP:0010850" %in% hpo_pat2$HPO
"HP:0002353" %in% hpo_pat1$HPO

mica("HP:0010850","HP:0002353")
ic_matrix["HP:0010850","HP:0002353"]
subset(allHPOs,allHPOs$term == "HP:0010850")$definition
subset(allHPOs,allHPOs$term == "HP:0002353")$definition

#both meet at EEG abnormality 

### -> ic_matrix is correct 

max_col <- apply(ic_matrix,2,max)

max_col %>% length() # 22 

max_row <- apply(ic_matrix,1,max)

max_row %>% length() # 27

max_complete <- sum(max_col,max_row)/2

#max_complete = 43.20424

round(max_complete,5) == 43.20424

