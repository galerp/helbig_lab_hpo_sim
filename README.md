# helbig_lab
This is a repository for collaborative work within the Helbig Lab research team and outside collaborators.


## Scripts


### IH_exact_p_HPO_joint_script_cleaned-2018-11-01.R
#### -Input: 
#####        * variant file - V7_Cumulative_328_Freeze_csv.csv
#####        * sim score


### MICA_Resnik_average_clean.R
#### -Input: 
#####       * Local IC - Local_IC_320_.csv
#####       * full HPO table - HPO_is.a_tree.csv
#####       * base table - Local_Base_per_patient_320.csv
#####       * prop table - Local_Prop_per_patient_320.csv
        
#### -Output: 
#####       * sim matrix
  

### MICA_Resnik_default_clean.R
#### -Input: 
#####       * Local IC - Local_IC_320_.csv
#####       * full HPO table - HPO_is.a_tree.csv
#####       * base table - Local_Base_per_patient_320.csv
#####       * prop table - Local_Prop_per_patient_320.csv      

###  calc_local_IC_2018-10-26.R
#### -Input: 
#####       * expanded file - 320_Expanded_v1_PG_clean.csv
#####       * full HPO table - HPO_is.a_tree.csv
#####       * HPO ancestors - HPO_Ancestors_Compressed.csv
#####       * n1_100k - 100k_resnik_def_log10_2018-10-30.csv
n1_100k
#### -Output: 
#####       * Propagated IC
#####       * Local IC - the local IC based on the propagated HPO term frequency in the cohort
#####       * Base Table - a stacked file of all the base terms from the cohort
#####       * Prop Table - a sacked file of all the propagated terms from the cohort
