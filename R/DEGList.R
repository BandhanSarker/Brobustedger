######################################################################################################
#######################################################################################################
#' @title Sample RNA-seq count data contains outlier values
#' @description sample RNAseq count \code{data-matrix} contains outlier values which is used for testing the package.
#' @format RNAseq count \code{data_matrix},whose row contains 1000 genes and column contains 6 samples, where 3 samples are control groups and 3 sample are case groups
#' @examples
#' data(sampleDataOut)
"sampleDataOut"

#' @title Sample RNA-seq count data contains missing values
#' @description sample RNAseq count \code{data-matrix} contains missing values which is used for testing the package.
#' @format RNAseq count \code{data_matrix},whose row contains 1000 genes and column contains 6 samples, where 3 samples are control groups and 3 sample are case groups
#' @examples
#' data(sampleDataMiss)
"sampleDataMiss"

##################################################################################
##################################################################################
######                           Boosting Robust edgeR                 ###########
##################################################################################
#' Identify differentially expressed gene list
#' @description DEGList is used to identify differentially expressed genes from RNA-seq count data.
#' @param data  RNA-seq count data matrix, whose row contain genes and column contain samples
#' @param n1    Control group
#' @param n2    Case group
#' @param p.threshold adjusted pvalue by default=0.05
#' @param padjust select the p-value adjustment method by default="BH". choose the possible method "none", "BH", "fdr", "BY" and "holm". See \code{\link[stats]{p.adjust}} for more details.
#' @param mi maximum number of iterations to impute nullify values by default= 5
#' @param NT number of trees in each set to impute nullify values by default=100
#' @return 1. Gene.list: Differentially expressed genes list and 2. result: others output
#' @import edgeR
#' @import limma
#' @import missRanger
#' @import stats
#' @importFrom methods is new
#' @seealso Outliers identification function: \code{\link[Brobustedger]{iLOO}}, imputaion of missing/nullify function: \code{\link[Brobustedger]{NI}}, See details \code{\link[edgeR]{DGEList}}, \code{\link[missRanger]{missRanger}}, \code{\link[Brobustedger]{sampleDataMiss}} and \code{\link[Brobustedger]{sampleDataOut}}
#' @examples
#' Brobustedger::DEGList(sampleDataOut,3,3,0.05) #Input RNA-seq count data
#' Brobustedger::DEGList(sampleDataMiss,3,3,0.05) #Input RNA-seq count data with missing or null value
#' @author Bandhan Sarker, Md. Matiur Rahaman, Muhammad Habibulla Alamin, Md. Ariful Islam, and Md. Nurul Haque Mollah
#' @references
#'Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#'Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
#'
#'George NI, Bowyer JF, Crabtree NM, Chang C-W (2015) An Iterative Leave-One-Out Approach to Outlier Detection in RNA-Seq Data. PLoS ONE 10(6): e0125224. https://doi.org/10.1371/journal.pone.0125224
#' @export
DEGList<-function(data,n1,n2,p.threshold=0.05,padjust="BH",mi=5,NT=100){

  ##############################
  message("Brobustedger::INFO:Computation takes few times... ...")
  #check_missing_values or not
  set.seed(100)
  if (sum(is.na(data)) > 0) {
    nmdata<- missRanger::missRanger(data, pmm.k = 3, maxiter = mi,
                                        num.trees = NT, verbose = 0)
   } else {
    nmdata<-data
  }
  #############
  sub<-nmdata

###############
  outData<-iLOO(sub)
  imputeData<-nmdata
  ######Test each gene whether it is an outlying gene or not using iLOO######
  location<-which(!is.na(outData))
  if(length(location)==0){
    dataM<-nmdata
  }else{
    ###Outlier’s value in each row and column should be considered as NULL#########

    imputeData[!is.na(outData)]<-NA
    dataM<-imputeData
  }


  ########Imputation of nullify values using RF############
  m<-NI(dataM)
  ## Added
  decideTestsDGE <- function(object,adjust.method="BH",p.value=0.05,lfc=0)
    #	Accept or reject hypothesis tests across genes and contrasts
    #	edgeR team. Original author was Davis McCarthy.
    #	Created 15 August 2010. Last modified 15 July 2018.
  {
    #	Check object class
    if( !(is(object,"DGEExact") || is(object,"DGELRT")) ) stop("Need DGEExact or DGELRT object")

    #	Apply multiple testing
    p <- object$table$PValue
    p <- p.adjust(p, method=adjust.method)
    isDE <- as.integer(p < p.value)

    #	Extract logFC
    logFC <- object$table$logFC

    #	Check for F-test with multiple logFC columns
    FTest <- is.null(logFC)

    #	With multiple contrasts, apply lfc threshold to maximum logFC
    if(FTest) {
      if(lfc>0) {
        coef.col <- grep("^logFC",colnames(object$table))
        logFC <- object$table[,coef.col]
        SmallFC <- rowSums(abs(logFC) >= lfc) == 0
        isDE[SmallFC] <- 0L
      }

      #	With single contrast, apply directionality and lfc threshold
    } else {
      isDE[isDE & logFC<0] <- -1L
      SmallFC <- (abs(logFC) < lfc)
      isDE[SmallFC] <- 0L
    }

    #	Assemble TestResults object
    isDE <- matrix(isDE, ncol=1)
    row.names(isDE) <- row.names(object)
    colnames(isDE) <- paste(rev(object$comparison),collapse="-")

    #	Record possible values
    if(FTest) {
      attr(isDE,"levels") <- c(0L,1L)
      attr(isDE,"labels") <- c("NotSig","Sig")
    } else {
      attr(isDE,"levels") <- c(-1L,0L,1L)
      attr(isDE,"labels") <- c("Down","NotSig","Up")
    }

    new("TestResults", isDE)
  }
  ########
  dataR<-round(m,0)
  # Assign condition
  control<-"control";case<-"case"
  condition <- factor(c(rep(control, n1), rep(case, n2)))
  coldata <- data.frame(row.names=colnames(dataR), condition)
  ################

  dge <- edgeR::DGEList(counts=dataR, group=coldata$condition)
  # Create the contrast matrix
  design.mat <-model.matrix(~ 0 + dge$samples$group)
  colnames(design.mat) <- levels(dge$samples$group)

  # Estimate dispersion parameter for GLM
  dge<-edgeR::estimateGLMRobustDisp(dge,  design.mat,maxit = 5, residual.type = "pearson")

  # Design matrix
  design.mat <-model.matrix(~ 0 + dge$samples$group)
  colnames(design.mat) <- c(control, case)
  # Model fitting

  fit.edgeR <- edgeR::glmFit(dge, design.mat)
  # Differential expression
  contrasts.edgeR <- limma::makeContrasts(control - case, levels=design.mat)
  lrt.edgeR <- edgeR::glmLRT(fit.edgeR, contrast=contrasts.edgeR)
  #logfold
  edgeR_results <- lrt.edgeR$table
  sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method=padjust, p.value = p.threshold)
  genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
  result<-edgeR_results[genes.edgeR,]
  result$p.adjust<-p.adjust(result$PValue, method =padjust, n = length(genes.edgeR))
  list(Gene.list=genes.edgeR,result=result)
}
#######################################################################################
################################## Outliers identification
#######################################################################################
#'@title Identify Outlier values from specific genes
#'@description iLOO is used to identify gene specific outliers from RNA-seq count data.
#'@param data  RNA-seq count data matrix, whose row contain genes and column contain samples
#'@return Output of the function identify gene specific outliers. Non-outlier values are displayed as "NA"
#'@seealso Differentially expressed gene list: \code{\link[Brobustedger]{DEGList}}, imputaion of missing/nullify function: \code{\link[Brobustedger]{NI}}, See details \code{\link[Brobustedger]{sampleDataMiss}} and \code{\link[Brobustedger]{sampleDataOut}}
#'@examples
#'Brobustedger::iLOO(sampleDataOut) # Input RNA-seq count data
#######
#' @export
iLOO <- function(data){
#estimate sequencing depth and compute cutoff
sub=data
sd <- mean(colSums(sub,na.rm=T))
sdcut <- 1/sd
suppressWarnings(
  #implement iterative scheme
  outl <- apply(sub, 1, function(z) {
    x <- z
    nbp <- 0
    out <- rep(0,0)
    track <- c(1:length(x))

    #iterative scheme
    while((min(nbp,na.rm=T) < sdcut) & (length(x)>2)) {

      #build matrix with rows representing leave-one out observation
      mat <- matrix(rep(x,length(x)),ncol=length(x),byrow=T)
      diag(mat) <- NA
      tmp <- t(mat)
      mat <- t(matrix(tmp[!is.na(tmp)],nrow=(length(x)-1),ncol=(length(x))))

      #fit negative binomial or Poisson distribution
      nbfit <- apply(mat,1,function(y) {
        if(length(y)>1) {
          v <- var(y)
          m <- mean(y)
          if(all(y==0)=="TRUE") {
            output <- NA
          } else if (v>m) {
            p <- mean(y)/var(y)
            r <- mean(y)^2/(var(y)-mean(y))
            output <- c(p,r)
          } else {
            lamb <- mean(y)
            output <- c(lamb)
          }
        } else output <- NA

        list(output)
      })

      nbfit <- lapply(nbfit, "[[", 1)
      #compute probabilities for leave-one out observation
      nbp <- rep(0,0)
      for (i in 1:length(nbfit)) {
        if(length(nbfit[[i]])==2) {
          nbp <- c(nbp,dnbinom(x[i],prob=nbfit[[i]][1],size=nbfit[[i]][2]))
        } else {
          nbp <- c(nbp,dpois(x[i],lambda=nbfit[[i]][1]))
        }
      }

      #compare probabilities to cutoff
      sel <- which(nbp < sdcut)
      if(length(sel)>0) x <- x[-sel]
      if(length(out)==0) {
        out <- c(out,track[sel])
      } else {
        out <- c(out,track[-out][sel])
      }

      fout <- rep(NA,length(z))
      if(length(out)>0) {
        fout[out] <- z[out]
      }
    }

    list(fout)
  })
)
#new data matrix with outliers (all other data is NA’ed)
identout <- matrix(unlist(outl),nrow=nrow(sub),ncol=ncol(sub),byrow=T)
colnames(identout) <- colnames(sub)
rownames(identout) <- rownames(sub)
return(identout)
}

#######################
#######################################################################################
##################################  identification of Null/missing value
#######################################################################################
#'@title NULL value Imputation
#'@description NI is used to impute null/missing values in the RNA-seq count data.
#'@param data  RNA-seq count data matrix contains missing/Null value, whose row contain genes and column contain samples
#'@param mi maximum number of iterations to impute nullify values
#'@param NT number of trees in each set to impute nullify values
#'@return imputation of missing/nullify data
#'@seealso Differentially expressed gene list: \code{\link[Brobustedger]{DEGList}}, Outliers identification function: \code{\link[Brobustedger]{iLOO}}, See details \code{\link[missRanger]{missRanger}} , \code{\link[Brobustedger]{sampleDataMiss}} and \code{\link[Brobustedger]{sampleDataOut}}
#'@examples
#'Brobustedger::NI(sampleDataMiss)  #Input RNA-seq count data with missing or null value
#' @export
NI <- function(data, mi = 5, NT = 100) {
  imputed <- missRanger::missRanger(data, pmm.k = 3, maxiter = mi,
                                    num.trees = NT, verbose = 0)
  return(imputed)
}
#######################################
#'_PACKAGE package
#'@name Brobustedger
#'@title Boosting Robust EdgeR
#'@keywords package
#'@description Brobustedger package is used to identify differentially expressed genes (DEGs) from RNA-Seq count data if the dataset contain outliars or missing values. It uses iLOO method for searching the location of outliers and then considered as missing values. Missing values are imputed using random forest method. Finally DEGs are calculated by robust edgeR.
#'@seealso Differentially expressed gene list: \code{\link[Brobustedger]{DEGList}}, Outliers identification function: \code{\link[Brobustedger]{iLOO}}, imputaion of missing/nullify function: \code{\link[Brobustedger]{NI}}, See details \code{\link[edgeR]{DGEList}} and \code{\link[missRanger]{missRanger}}
#'@author Bandhan Sarker, Md. Matiur Rahaman, Muhammad Habibulla Alamin, Md. Ariful Islam, and Md. Nurul Haque Mollah
#'@references
#'Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#'Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
#'
#'George NI, Bowyer JF, Crabtree NM, Chang C-W (2015) An Iterative Leave-One-Out Approach to Outlier Detection in RNA-Seq Data. PLoS ONE 10(6): e0125224. https://doi.org/10.1371/journal.pone.0125224
NULL
#####################################################################################
####                                END                                   ###########
#####################################################################################
