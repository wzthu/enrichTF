#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlapPairs
#' @importFrom S4Vectors second first
#' @importFrom utils write.table read.table
#'
setClass(Class = "RegionConnectTargetGene",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        inputForegroundBed <- allparam[["inputForegroundBed"]]
        inputBackgroundBed <- allparam[["inputBackgroundBed"]]
        outputForegroundBed <- allparam[["outputForegroundBed"]]
        outputBackgroundBed <- allparam[["outputBackgroundBed"]]
        regularGeneCorrBed <- allparam[["regularGeneCorrBed"]]
        enhancerRegularGeneCorrBed <- allparam[["enhancerRegularGeneCorrBed"]]
        ouputForgroundGeneTxt <- allparam[["ouputForgroundGeneTxt"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            input(.Object)$inputForegroundBed <- output(prevStep)$outputForegroundBed
            input(.Object)$inputBackgroundBed <- output(prevStep)$outputBackgroundBed
        }
        if(!is.null(inputForegroundBed)){
            input(.Object)$inputForegroundBed <- inputForegroundBed
        }
        if(!is.null(inputBackgroundBed)){
            input(.Object)$inputBackgroundBed <- inputBackgroundBed
        }



        if(is.null(outputForegroundBed)){
            output(.Object)$outputForegroundBed <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputForegroundBed"]],
                            regexSuffixName = "foreground.bed",
                            suffix = "gene.foreground.bed")
        }else{
            output(.Object)$outputForegroundBed <- outputForegroundBed
        }

        if(is.null(ouputForgroundGeneTxt)){
            output(.Object)$ouputForgroundGeneTxt <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputForegroundBed"]],
                            regexSuffixName = "foreground.bed",
                            suffix = "gene.foreground.txt")
        }else{
            output(.Object)$ouputForgroundGeneTxt <- ouputForgroundGeneTxt
        }


        if(is.null(outputBackgroundBed)){
            output(.Object)$outputBackgroundBed <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputBackgroundBed"]],
                            regexSuffixName = "background.bed",
                            suffix = "gene.background.bed")
        }else{
            output(.Object)$outputBackgroundBed <- outputBackgroundBed
        }

        if(is.null(regularGeneCorrBed)){
            output(.Object)$regularGeneCorrBed <-
                getRefFiles("RE_gene_corr")
        }else{
            output(.Object)$regularGeneCorrBed <-
                regularGeneCorrBed
        }

        if(is.null(enhancerRegularGeneCorrBed)){
            output(.Object)$enhancerRegularGeneCorrBed <-
                getRefFiles("Enhancer_RE_gene_corr")
        }else{
            output(.Object)$enhancerRegularGeneCorrBed <-
                enhancerRegularGeneCorrBed
        }
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,...){
        inputForegroundBed <- getParam(.Object,"inputForegroundBed")
        inputBackgroundBed <- getParam(.Object,"inputBackgroundBed")
        outputForegroundBed <- getParam(.Object,"outputForegroundBed")
        outputBackgroundBed <- getParam(.Object,"outputBackgroundBed")
        regularGeneCorrBed <- getParam(.Object,"regularGeneCorrBed")
        enhancerRegularGeneCorrBed <- getParam(.Object,
                                               "enhancerRegularGeneCorrBed")
        ouputForgroundGeneTxt <- getParam(.Object, "ouputForgroundGeneTxt")

        inputForegroundgr <- import(con=inputForegroundBed)
        inputBackgroundgr <- import(con=inputBackgroundBed)

        rg <- import(con=regularGeneCorrBed,
                     colnames=c("name","score","blockCount"))
        erg <- import(con=enhancerRegularGeneCorrBed,
                      colnames=c("name","score","blockCount"))

        pairs <- findOverlapPairs(inputForegroundgr,rg)
        first(pairs)$geneName  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputForegroundgr <- first(pairs)

        pairs <- findOverlapPairs(inputForegroundgr,erg)
        first(pairs)$geneName  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputForegroundgr <- c(outputForegroundgr,first(pairs))

        pairs <- findOverlapPairs(inputBackgroundgr,rg)
        first(pairs)$geneName  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputBackgroundgr <- first(pairs)

        pairs <- findOverlapPairs(inputBackgroundgr,erg)
        first(pairs)$geneName  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputBackgroundgr <- c(outputBackgroundgr,first(pairs))

        write.table(as.data.frame(outputForegroundgr[
            mcols(outputForegroundgr)$score>0.3])
            [,c("seqnames","start","end","name",
                "score","geneName","blockCount")],
            outputForegroundBed,sep="\t",quote = FALSE,
            row.names = FALSE,col.names = FALSE)

        write.table(as.data.frame(geneName = mcols(outputForegroundgr[
            mcols(outputForegroundgr)$score>0.3])$geneName
            ),
            ouputForgroundGeneTxt,sep="\t",quote = FALSE,
            row.names = FALSE,col.names = FALSE)

        write.table(as.data.frame(outputBackgroundgr[
            mcols(outputBackgroundgr)$score>0.3])
            [,c("seqnames","start","end","name",
                "score","geneName","blockCount")],
                    outputBackgroundBed,sep="\t",quote = FALSE,
            row.names = FALSE,col.names = FALSE)


#        export.bed(outputForegroundgr,outputForegroundBed)
#        export.bed(outputBackgroundgr,outputBackgroundBed)


#        .Object@propList[["motifs_in_region"]] <- result

        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,...){
        if(is.null(.Object$inputForegroundBed)){
            stop("inputForegroundBed is required.")
        }
        if(is.null(.Object$inputBackgroundBed)){
            stop("inputBackgroundBed is required.")
        }

    }
)



setMethod(
    f = "checkAllPath",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,...){
        checkFileExist(.Object$inputForegroundBed)
        checkFileExist(.Object$inputBackgroundBed)

    }
)


#' @name RegionConnectTargetGene
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
#' @param inputForegroundBed \code{Character} scalar.
#' The BED file directory of foreground regions.
#' @param inputBackgroundBed  \code{Character} scalar.
#' The BED file directory of background regions.
#' @param outputForegroundBed \code{Character} scalar.
#' The BED file directory of target genes connecting with foreground regions,
#' which are derived from PECA model.
#' Default: NULL (generated base on inputForegroundBed)
#' @param outputBackgroundBed \code{Character} scalar.
#' The BED file directory of target genes connecting with background regions,
#' which are derived from PECA model.
#' Default: NULL (generated base on inputBackgroundBed)
#' @param regularGeneCorrBed \code{Character} scalar.
#' The BED file directory of target genes which are predicted from PECA.
#' Default: NULL (e.g. after \code{library (enrichTF)}, you can call
#' function \code{setGenome("hg19")})
#' @param enhancerRegularGeneCorrBed \code{Character} scalar.
#' The BED file directory of enhancer-targets predicted from PECA.
#' Default: NULL (e.g. after \code{library (enrichTF)}, you can call
#' function \code{setGenome("hg19")})
#' @param ouputForgroundGeneTxt \code{Character} scalar.
#' The TXT file directory of target genes list connecting with foreground regions,
#' which are derived from PECA model.
#' Default: NULL (generated base on inputForegroundBed)
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
#' conTG <- enrichRegionConnectTargetGene(gen)



setGeneric("enrichRegionConnectTargetGene",function(prevStep,
                                                    inputForegroundBed = NULL,
                                                    inputBackgroundBed = NULL,
                                                    outputForegroundBed = NULL,
                                                    outputBackgroundBed = NULL,
                                                    regularGeneCorrBed = NULL,
                                                    enhancerRegularGeneCorrBed = NULL,
                                                    ouputForgroundGeneTxt = NULL,
                                                    ...) standardGeneric("enrichRegionConnectTargetGene"))



#' @rdname RegionConnectTargetGene
#' @aliases enrichRegionConnectTargetGene
#' @export
setMethod(
    f = "enrichRegionConnectTargetGene",
    signature = "Step",
    definition = function(prevStep,
                          inputForegroundBed = NULL,
                          inputBackgroundBed = NULL,
                          outputForegroundBed = NULL,
                          outputBackgroundBed = NULL,
                          regularGeneCorrBed = NULL,
                          enhancerRegularGeneCorrBed = NULL,
                          ouputForgroundGeneTxt = NULL,
                          ...){
        allpara <- c(list(Class = "RegionConnectTargetGene",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname RegionConnectTargetGene
#' @aliases regionConnectTargetGene
#' @export
regionConnectTargetGene <- function(inputForegroundBed,
                                    inputBackgroundBed,
                                    outputForegroundBed = NULL,
                                    outputBackgroundBed = NULL,
                                    regularGeneCorrBed = NULL,
                                    enhancerRegularGeneCorrBed = NULL,
                                    ouputForgroundGeneTxt = NULL,
                                    ...){
    allpara <- c(list(Class = "RegionConnectTargetGene", prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
