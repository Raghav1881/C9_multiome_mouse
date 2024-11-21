#  Integrative Analysis of Single-Nucleus RNA-Seq and ATAC-Seq of C9orf72 Knockout Mouse Brain
![alt txt](https://github.com/user-attachments/assets/6b216559-61d7-4f75-99d1-be5d9df5c233)

This repository contains code for the analysis of C9orf72 haploinsufficiency and knockout in the hippocampus and frontal cortex.

## Abstract
Hexanucleotide repeat expansions in the first intron of C9orf72 are the most frequent genetic cause of amyotrophic lateral sclerosis (ALS) and frontotemporal dementia (FTD). These expansions result in neurodegeneration related to haploinsufficiency of C9orf72, though the mechanisms remain unclear. To investigate the effects of C9orf72 haploinsufficiency, I conducted single nucleus RNA-seq (snRNA-seq) and single nucleus assay for transposase-accessible chromatin sequencing (snATAC-seq) on hippocampus and frontal cortex tissues from C9orf72 heterozygous-knockout (C9-HET) and homozygous-knockout (C9-KO) mice. In the hippocampus, C9-HET mice showed downregulation of excitotoxicity-related pathways in neurons, whereas C9-KO mice showed upregulation. In the frontal cortex, C9-KO showed dysregulation of protein folding, ATP synthesis, and epigenomic dysregulation in glial cells and layer 2-3 excitatory neurons, highlighting a potential vulnerability of specific cell types. This multiome analysis provides insights into the molecular mechanisms underlying C9orf72-related ALS and FTD.

## Folder Structure 
There are two main folders: C9Mouse_Frontal_Cortex and hippocampus_multiome_2024 containing analyses of the frontal cortex and hippocampus respectively. Within each, there is a main.R file which sources all functions used from the R folder to run analyses. Other folders correspond to outputs and objects/figures created.

In the parent folder, there are two jupyter notebooks: ccans.ipnyb and figure_plot.ipnyb which contain code to generate ccans and figures used in the thesis respectively. There is also a requirements.txt that can be used to create a conda environment.  

## Data Availability
[Thesis can be found here.](https://hdl.handle.net/1807/141248)
Fastq files are available from the SRA with the following bioproject: PRJNA1166670 

Seurat objects for the hippocampus and frontal cortex can be made available upon request.

## References
If any data or code is used from this analysis, please cite using the following:

Sharma, R. Integrative Analysis of Single-Nucleus RNA-Seq and ATAC-Seq of C9orf72 Knockout Mouse Brain: Unraveling the Molecular Mechanisms Underlying Neurodegeneration in Amyotrophic Lateral Sclerosis and Frontotemporal Dementia. (2024).

