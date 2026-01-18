# Brobustedger
**Brobustedger-package** is developed to identify differentially expressed genes (DEGs) from RNA-seq count data through **boosting robust edgeR**  when the dataset contains outlier/noise or missing values. Iterative leave-out (iLOO) method is used to detect of gene specific outliers and random forest is used to impute null/missing values. Finally DEGs are calculated by edgeR (Robust). **Brobustedger-package** provides a user-friendly R-package interface.  Users can utilize its **DEGList** function to find differentially expressed gene list. Additionally, **Brobustedger-package** offers functionalities for outlier detection (**iLOO**) and imputation of missing or null values (**NI**) by using *iLOO* and *random forest* method, respectively.   

## Updates
- The package now uses the `missRanger` package instead of `missForest` for missing and outlier value imputation, providing faster performance and better
  scalability for large datasets.

## Installation
```r
#Install edgeR from bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
# Install devtools from CRAN
install.packages("devtools")

# development version from GitHub:
require("devtools")
devtools::install_github("BandhanSarker/Brobustedger",build_vignettes = TRUE)
```
## Data Pre-process
## Step-1:
iLOO is used to identify gene specific outliers from RNA-seq count data.
## Usage
```r
iLOO(data)
```
## Arguments
* **data** RNA-seq count data matrix, whose row contain genes and column contain samples
## Value
Output of the function identify gene specific outliers. Non-outlier values are displayed as "NA"
## Examples
```r
Brobustedger::iLOO(sampleDataOut) # Input RNA-seq count data
```
## Step-2:
NI is used for the imputation of null/missing values in RNA-seq count data.
## Usage
```r
NI(data, mi = 5, NT = 100)
```
## Arguments
* **data**	RNA-seq count data matrix contains missing/Null value, whose row contain genes and column contain samples
* **mi**	maximum number of iterations to impute nullify values
* **NT**    number of trees in each set to impute nullify values

## Value
Clean dataset after the imputation of missing/nullify value

## Examples
```r
Brobustedger::NI(sampleDataMiss)  #Input RNA-seq count data with missing or null value
```
## Step-3:
DEGList is used to identify differentially expressed genes from RNA-seq count data.
## Usage
```r
#Identify differentially expressed gene list
DEGList (data, n1, n2, p.threshold = 0.05)
```
## Arguments
* **data**	  RNA-seq count data matrix, whose row contains genes and column contains  samples                 
* **n1**	  control group
* **n2**	  case group
* **p.threshold**	adjusted p-value by default=0.05
* **padjust**	select the p-value adjustment method by default="BH". choose the possible method "none", "BH", "fdr", "BY" and "holm". See help section of BrobustedgeR-package for more details.
* **mi** maximum number of iterations to impute nullify values by default= 5
* **NT** number of trees in each set to impute nullify values by default=10
  
## Value 
1. Gene.list: Differentially expressed genes list and 2. result: others output
## Examples
```r
 data (sampleDataOut)
 DEGList(sampleDataOut,3,3,0.05)
```
## Example Dataset
```r
#Sample RNA-seq count data with outliers value
data(sampleDataOut)
```
```r
#Sample RNA-seq count data with missing value
data(sampleDataMiss)
```
 
## Author(s): 
Bandhan Sarker, Md. Matiur Rahaman, Muhammad Habibulla Alamin, Md. Ariful Islam and Md.  Nurul Haque Mollah

## References
Sarker B, Rahaman MM, Alamin MH et al. (2024) Boosting edgeR (Robust) by dealing with missing observations and gene-specific outliers in RNA-Seq profiles and its application to explore biomarker genes for diagnosis and therapies of ovarian cancer. Genomics. https://doi.org/10.1016/j.ygeno.2024.110834

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140

Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597

George NI, Bowyer JF, Crabtree NM, Chang C-W (2015) An Iterative Leave-One-Out Approach to Outlier Detection in RNA-Seq Data. PLoS ONE 10(6): e0125224. https://doi.org/10.1371/journal.pone.0125224

