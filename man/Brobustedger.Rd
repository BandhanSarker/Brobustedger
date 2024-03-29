% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEGList.R
\docType{package}
\name{Brobustedger}
\alias{Brobustedger-package}
\alias{Brobustedger}
\title{Boosting Robust EdgeR}
\description{
Brobustedger package is used to identify differentially expressed genes (DEGs) from RNA-Seq count data if the dataset contain outliars or missing values. It uses iLOO method for searching the location of outliers and then considered as missing values. Missing values are imputed using random forest method. Finally DEGs are calculated by robust edgeR.
}
\references{
Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140

Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597

George NI, Bowyer JF, Crabtree NM, Chang C-W (2015) An Iterative Leave-One-Out Approach to Outlier Detection in RNA-Seq Data. PLoS ONE 10(6): e0125224. https://doi.org/10.1371/journal.pone.0125224
}
\seealso{
Differentially expressed gene list: \code{\link[Brobustedger]{DEGList}}, Outliers identification function: \code{\link[Brobustedger]{iLOO}}, imputaion of missing/nullify function: \code{\link[Brobustedger]{NI}}, See details \code{\link[edgeR]{DGEList}} and \code{\link[missForest]{missForest}}
}
\author{
Bandhan Sarker, Md. Matiur Rahaman, Muhammad Habibulla Alamin, Md. Ariful Islam, and Md. Nurul Haque Mollah
}
\keyword{package}
