# HPO_Sim_Gene

Using the Human Phenotype Ontology (HPO) and a cohort of patients annotated with HPO terms and VCF files, these scripts find phenotype-genotype correlations. This can aid in gene discovery, treatment, and a better understanding of genes' phenotypic variability.
First, using HPO terms similarity scores are found for every patient pair. Next, genes with potentially causitive variants in multiple patients are extracted from the VCF and the median similarity scores among each of these patients is calculated.  Using permutation analysis of median similarity scores, p-values are assigned to each of these genes. A lower p-value potentially indicates a causitive gene.

 ### Requirements:
  
 * [R](https://www.r-project.org/) with packages tidyverse and memoise.
                
 
 
 ### Steps to Run:

 * Clone the repository, modify the [config file](https://github.com/galerp/hpo_sim_gene/blob/master/input.yml).
 
 * In the [config file](https://github.com/galerp/hpo_sim_gene/blob/master/input.yml) mention the the field `output_dir`, this is where your output files would be written to. 

 * Run [R file](https://github.com/galerp/hpo_sim_gene/blob/master/hpo_dee_sim.R), specifying the YAML config file using the `--input` flag .
```
~/Rscript hpo_dee_sim.R --input /path_to/input.yml
```
 
 
 ### Running the tests

There are test files available here: [Files](https://github.com/galerp/hpo_sim_gene/tree/master/Files). 
Ensure that these files are linked appropriately in the [config file](https://github.com/galerp/hpo_sim_gene/blob/master/input.yml) as such:
```
patient_phenome : Files/example_phenome.csv  

variant_file : Files/variant.csv  
```
This provides the necessary VCF and cohort of patients with annotated HPO terms.
