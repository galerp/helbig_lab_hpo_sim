#compare atonic seizures and spasms

#HPO atonic seizures 
# look up atonic seizures in allHPOs
# HP:0010819
# check whether this is present in the local_IC
# go to local_IC
# yes, row 7489

#HPO infantile spasms
# look up infantile spasms in allHPOs
# HP:0012469
# check whether this is present in the local_IC
# go to local_IC
# yes, row 9121

#run mica function
overlap = mica("HP:0010819","HP:0012469")

hpo1 = "HP:0010819" # atonic seizures
hpo2 = "HP:0012469" # infantile spasms

#taken from MICA function
path1_unique <- path[which(path$Term == hpo1),]

path2_unique <- path[which(path$Term == hpo2),]

#add term definition
path1_unique %>% left_join(allHPOs, by=c('HPO1'='term')) -> path1_def 
path2_unique %>% left_join(allHPOs, by=c('HPO1'='term')) -> path2_def 

#look at all common ancestors
joint1 <- path1_unique %>% inner_join(path2_unique, by = 'Counter1') 

#add IC to find MICA 
joint2 <- joint1 %>% left_join(ic2, by = c('Counter1' = 'HPO')) 

#get MICA value
mica_ic <- joint2$Propagated.local.IC %>% max  

#should be MICA for "seizures" - HP:0001250 
#go to local_IC for HP:0001250
#propagated local IC = 0.01120822
#mica_ic = 0.01120822
mica_ic


