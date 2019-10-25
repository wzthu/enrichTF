
setClass(Class = "TissueOpennessSpecificity",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "TissueOpennessSpecificity",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        bedOutput <- allparam[["bedOutput"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            bedInput <- getParam(prevStep,"bedOutput")
            input(.Object)$bedInput <- foregroundBed
        }

        if(!is.null(bedInput)){
            input(.Object)$bedInput <- bedInput
        }

        if(is.null(bedOutput)){
            output(.Object)$bedOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed|bed.gz|bed.bz2",suffix = "bed")
        }else{
            output(.Object)$bedOutput <- bedOutput
        }

        .Object
    }
)



setMethod(
    f = "processing",
    signature = "UnzipAndMergeBed",
    definition = function(.Object,...){
        bedInput <- getParam(.Object,"bedInput")
        bedOutput <- getParam(.Object,"bedOutput")

        beds <- lapply(bedInput, function(bedfile){
            stopifnot(is.character(bedfile))
            stopifnot(file.exists(bedfile))
            tmpfile <- getAutoPath(.Object, bedfile, "bed|bed.gz|bed.bz2","tmp.bed")
            decompressBed(.Object, bedfile, tmpfile)
            bed <- read.table(tmpfile, header = FALSE, sep = "\t")
            bed <- bed[,1:3]
            colnames(bed) <-c("chrom", "start", "end")
            unlink(tmpfile)
            return(bed)
        })

        bed <- do.call(c, beds)

        bed <- as(bed, "GRanges")

        export.bed(bed, con = bedOutput)

        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "UnzipAndMergeBed",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["bedInput"]])){
            stop("bedInput is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "UnzipAndMergeBed",
    definition = function(.Object,...){
        checkFileExist(input(.Object)[["bedInput"]]);

    }
)


#' @name UnzipAndMergeBed
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Generate background regions and reset the size of
#' foreground regions
#' @description
#' Use uniform distribution to generate background
#' sequence regions from genome.
#' The size of foreground regions will be unified into the
#' length specified in argument.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process
#' from other packages, such as esATAC.
#' @param bedInput \code{Character} scalar or vector.
#' The directory of region BED files for analysis.
#' BED, BED.gz, BED.bz2 formats are supported.
#' @param bedOutput \code{Character} scalar.
#' The BED output file directory of merged BED files.
#' Default: NULL (generated base on bedInput)
#' @param ... Additional arguments, currently unused.
#' @details
#' All compressed files will be de-compressed.
#' Only first 3 columns (chromasomes, start and end) will be collected.
#' All BED files will be merged into one BED file.
#' @return An invisible \code{\link{EnrichStep-class}}
#' object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{genBackground}}


#' @examples
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed.gz")
#' gen <- unzipAndMergeBed(bedInput = foregroundBedPath)


setGeneric("enrichUnzipAndMergeBed",
           function(prevStep,
                    bedInput = NULL,
                    bedOutput = NULL,
                    ...) standardGeneric("enrichUnzipAndMergeBed"))



#' @rdname TissueOpennessSpecificity
#' @aliases enrichTissueOpennessSpecificity
#' @export
setMethod(
    f = "enrichUnzipAndMergeBed",
    signature = "Step",
    definition = function(prevStep,bedInput, bedOutput = NULL, ...){
        allpara <- c(list(Class = "UnzipAndMergeBed",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname TissueOpennessSpecificity
#' @aliases tissueOpennessSpecificity
#' @export
unzipAndMergeBed <- function(bedInput, bedOutput = NULL, ...){
    allpara <- c(list(Class = "UnzipAndMergeBed",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
