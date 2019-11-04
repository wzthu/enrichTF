
setClass(Class = "GeneOntology",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "GeneOntology",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        inputTxt <- allparam[["inputTxt"]]
        outputTxt <- allparam[["outputTxt"]]
        outputPdf <- allparam[["outputPdf"]]
        orgDb <- allparam[["orgDb"]]
        keyType <- allparam[["keyType"]]

        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            input(.Object)$inputTxt <- output(prevStep)$ouputForgroundGeneTxt
        }
        if(!is.null(inputTxt)){
            input(.Object)$inputTxt <- inputTxt
        }

        if(!is.null(orgDb)){
            param(.Object)$orgDb <- orgDb
        }else{
            param(.Object)$orgDb <- getRefRc("OrgDb")
        }

        stopifnot(is.character(keyType))
        param(.Object)$keyType <- keyType


        if(is.null(outputTxt)){
            output(.Object)$outputTxt <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputTxt"]],
                            regexSuffixName = "txt",
                            suffix = "txt")
        }else{
            output(.Object)$outputTxt <- outputTxt
        }


        if(is.null(outputPdf)){
            output(.Object)$outputPdf <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputTxt"]],
                            regexSuffixName = "bed",
                            suffix = "pdf")
        }else{
            output(.Object)$outputPdf <- outputPdf
        }

        .Object
    }
)

#' @importFrom clusterProfiler enrichGO

setMethod(
    f = "processing",
    signature = "GeneOntology",
    definition = function(.Object,...){
        inputTxt <- getParam(.Object,"inputTxt")
        outputTxt <- getParam(.Object,"outputTxt")
        outputPdf <- getParam(.Object,"outputPdf")
        orgDb <- getParam(.Object,"orgDb")
        keyType <- getParam(.Object,"keyType")

        genelist <- read.table(inputTxt, header = FALSE, sep = '\t')
        genelist <- unique(as.character(genelist[,1]))

        ego<-clusterProfiler::enrichGO(gene = genelist,
                                       OrgDb = orgDb,
                                       keyType   = keyType)

        write.table(ego@result, file = outputTxt, sep='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
        file.create(outputPdf)
        .Object
    }
)



#' @name GeneOntology
#' @importFrom TFBSTools getMatrixSet
#' @importFrom TFBSTools PFMatrixList
#' @importFrom TFBSTools PFMatrix
#' @title Gene ontology enrichment analysis for provided gene list
#' @description
#' User provide target gene list of forground region.
#' This function will call function for gene ontology enrichment analysis
#' @param prevStep \code{\link{Step-class}} object scalar.
#' This parameter is available when the upstream step function
#' (printMap() to see the previous functions)
#' have been sucessfully called.
#' Accepted value can be the object return by any step function or be feed by
#' \code{\%>\%} from last step function.
#' @param inputTxt \code{Character} scalar.
#' Gene list text file. All gene names are in one column.
#' @param outputTxt  \code{Character} scalar.
#' Gene ontology enrichment analysis result table.
#' Each row contain one gene ontology information.
#' @param outputPdf \code{Character} scalar.
#' Gene ontology enrichment analysis result figure.
#' It contains gene ontology network.
#' @param orgDb \code{Character} scalar.
#' Bioconductor OrgDb object name for gene ontology enrichment analysis.
#' @param keyType \code{Character} scalar.
#' Gene name type include "SYMBOL" and "ENSEMBLE"
#' @param ... Additional arguments, currently unused.
#' @details
#' Currently, this function call enrichGO from package clusterProfiler to implement this funtion.
#' @return An invisible \code{\link{EnrichStep-class}} object
#' (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{regionConnectTargetGene}}
#' \code{\link{enrichRegionConnectTargetGene}}
#' @examples
#'
#' genelist.txt <- system.file(package = "enrichTF", "extdata","genelist.txt")
#' geneOntology(inputTxt = genelist.txt, orgDb = "org.Hs.eg.db")
#'



setGeneric("enrichGeneOntology",function(prevStep,
                                         inputTxt = NULL,
                                         outputTxt = NULL,
                                         outputPdf = NULL,
                                         orgDb = NULL,
                                         keyType = "SYMBOL",
                                         ...) standardGeneric("enrichGeneOntology"))



#' @rdname GeneOntology
#' @aliases enrichGeneOntology
#' @export
setMethod(
    f = "enrichGeneOntology",
    signature = "Step",
    definition = function(prevStep,
                          inputTxt = NULL,
                          outputTxt = NULL,
                          outputPdf = NULL,
                          orgDb = NULL,
                          keyType = "SYMBOL",
                          ...){
        allpara <- c(list(Class = "GeneOntology",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname GeneOntology
#' @aliases geneOntology
#' @export
geneOntology <- function(inputTxt,
                         outputTxt = NULL,
                         outputPdf = NULL,
                         orgDb = NULL,
                         keyType = "SYMBOL",
                         ...){
    allpara <- c(list(Class = "GeneOntology", prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
