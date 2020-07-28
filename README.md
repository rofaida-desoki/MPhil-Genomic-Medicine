# Scripts of the pipeline implemented in the MPhil research project

## MPhil in Genomic Medicine 2019-2020
## University of Cambridge

### Research Project Title:
Exploring Cancer Predisposition due to E-Cadherin-related and Mismatch Repair genes utilizing the 100,000 Genomes Project dataset

### Study Aims:

The identification of germline variants predisposing to cancer plays an essential role in exploring the molecular pathogenesis of tumours and is a critical step towards the clinical management of affected families.

The aim of this association analysis was to explore the 100K Genomes Project dataset for possible associations of germline variants in the selected cancer-associated pathways. 
The implemented methods includes Burden and Variance-componenet analysis of aggregated variants per gene and aggregated analysis of MMR pathway and HDGC-associated genes.
The 100K GP dataset provided a considerable amount of data and an equipped environment for the research purpose.

The selected cancer-associated pathways were as follows:
-	MMR pathway: The heterozygous pathogenic variants in MMR genes (MLH1, MSH2, MSH6, PMS2, MSH3, PMS1 and EPCAM) predispose to Lynch Syndrome. 
-	HDGC-related genes: Hereditary diffuse gastric cancer (HDGC) is an autosomal dominant cancer predisposition syndrome which is associated with pathogenic variants of CDH1, CTNNA1, MAP3K6, MYD88 genes.

## Study objectives: 

1. Explore the Genomics England (GEL) research environment.
2. Extract data from the latest main programme data release in LabKey server available in the Research Environment.
3. Select the sample cohort across multiple somatic cancers and cancer predisposition syndromes from the 100K GP participants.
4. Identify germline variants in the selected genes (MMR genes and HDGC-predisposing genes) across the selected cohorts.
5. Implement a range of statistical analysis techniques to analyse variants aggregated per gene and per whole pathway analysis.
6. Confirm previously established associations of genes and pathways with related types of cancer.
7. Examine if associations of each gene/pathway could be established with other types of cancer.

## The main steps in the pipeline are to:

1. R import: Import, explore and prepare participants data from LabKey
2. Sample selection: Select cases and controls
3. Genes selection: Identify gene coordinates
4. VCF preparation: Extract sequencing data from gVCF
5. Annotation: Ensembl VEP, CADD and ClinVar variants annotation 
6. QC and filtering: QC checks and functional filtering of variants
7. Data consolidation: VCF and phenotypic data into R and additioal genotype filtering
8. SKAT preparation: prepare files required for SKAT library 
9. Association analysis: Test for associations using SKAT library

The repository was set-up to support the MPhil research project through providing the details of the project pipeline.
The repository can be used by the public and it is adapted to the Genomics England Research Environment.

