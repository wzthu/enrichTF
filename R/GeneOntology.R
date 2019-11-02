
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
                            regexSuffixName = "bed",
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

setMethod(
    f = "checkRequireParam",
    signature = "GeneOntology",
    definition = function(.Object,...){

    }
)



setMethod(
    f = "checkAllPath",
    signature = "GeneOntology",
    definition = function(.Object,...){


    }
)


#' @name GeneOntology
#' @importFrom TFBSTools getMatrixSet
#' @importFrom TFBSTools PFMatrixList
#' @importFrom TFBSTools PFMatrix
#' @title Connect regions with their target genes
#' @description
#' Connect foreground and background regions to their target genes,
#' which is predicted from PECA model.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process from
#' \code{\link{genBackground}} or \code{\link{enrichGenBackground}}
#' when it is not used in a pipeline.  If it is used in a pipeline or
#' \code{\%>\%} is applied on this function, any steps in this package
#' is acceptable.
#' @param inputTxt \code{Character} scalar.
#' The BED file directory of foreground regions.
#' @param outputTxt  \code{Character} scalar.
#' The BED file directory of background regions.
#' @param outputPdf \code{Character} scalar.
#' The BED file directory of target genes connecting with foreground regions,
#' which are derived from PECA model.
#' Default: NULL (generated base on inputForegroundBed)
#' @param orgDb \code{Character} scalar.
#' The BED file directory of target genes connecting with background regions,
#' which are derived from PECA model.
#' Default: NULL (generated base on inputBackgroundBed)
#' @param keyType \code{Character} scalar.
#' The BED file directory of target genes which are predicted from PECA.
#' Default: NULL (e.g. after \code{library (enrichTF)}, you can call
#' function \code{setGenome("hg19")})
#' @param ... Additional arguments, currently unused.
#' @details
#' Connect foreground and background regions to target genes,
#' which are predicted from PECA.
#' @return An invisible \code{\link{EnrichStep-class}} object
#' (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{genBackground}}
#' \code{\link{findMotifsInRegions}}
#' \code{\link{tfsEnrichInRegions}}
#' @examples
#' setGenome("testgenome") #Use "hg19","hg38",etc. for your application
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' gen <- genBackground(inputForegroundBed = foregroundBedPath)
#' conTG <- enrichGeneOntology(gen)



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
