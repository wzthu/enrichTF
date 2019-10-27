
setClass(Class = "TissueOpennessConserve",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "TissueOpennessConserve",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        openConserveBedInput <- allparam[["openConserveBedInput"]]
        bedOutput <- allparam[["bedOutput"]]
        distrPdfOutput <- allparam[["distrPdfOutput"]]

        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            bedInput0 <- getParam(prevStep,"bedOutput")
            input(.Object)$bedInput <- bedInput0
        }

        if(!is.null(bedInput)){
            input(.Object)$bedInput <- bedInput
        }

        if(!is.null(openConserveBedInput)){
            input(.Object)$openConserveBedInput <- openConserveBedInput
        }else{
            input(.Object)$openConserveBedInput <- getRefFiles("ConserveRegion")
        }


        if(is.null(bedOutput)){
            output(.Object)$bedOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "bed")
        }else{
            output(.Object)$bedOutput <- bedOutput
        }

        if(is.null(distrPdfOutput)){
            output(.Object)$distrPdfOutput <-
                getAutoPath(.Object,originPath =
                                input(.Object)[["bedInput"]][1],
                            regexSuffixName = "bed",suffix = "pdf")
        }else{
            output(.Object)$distrPdfOutput <- distrPdfOutput
        }

        .Object
    }
)



setMethod(
    f = "processing",
    signature = "TissueOpennessConserve",
    definition = function(.Object,...){
        bedInput <- getParam(.Object,"bedInput")
        bedOutput <- getParam(.Object,"bedOutput")
        openConserveBedInput <-getParam(.Object,"openConserveBedInput")
        distrPdfOutput <- getParam(.Object,"distrPdfOutput")


        # read input data
        region <- import.bed(bedInput)
        openTable <- read.table(openConserveBedInput,header = F,sep = "\t")

        # transform open table into GRanges format
        openBed <- openTable[,1:3]
        colnames(openBed) <- c("chrom", "start", "end")
        openValue <- openTable[,4:ncol(openTable)]
        colnames(openValue) <- seq_len(ncol(openValue))
        openRanges <- as(openBed,"GRanges")
        mcols(openRanges) <- openValue

        # find overlapped region
        pairs <- findOverlapPairs(openRanges, region, ignore.strand = TRUE)
        openRegion <- first(pairs)
        openRegion <- as.data.frame(openRegion)
        write.table(openRegion[,c(1:3,6:ncol(openRegion))],file = bedOutput,
                    sep="\t", col.names = FALSE, row.names = FALSE)


        # draw distribution


        pdf(distrPdfOutput)
        ggplot(openRegion, aes(X2)) + geom_histogram(binwidth = 0.01) + xlab("conserve") + ylab("count")
        dev.off()


        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "TissueOpennessConserve",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["bedInput"]])){
            stop("bedInput is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "TissueOpennessConserve",
    definition = function(.Object,...){
        checkFileExist(input(.Object)[["bedInput"]]);

    }
)


#' @name TissueOpennessConserve
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Tissue's open conservation of the given region
#' @description
#' User provide region through a BED file.
#' This function will provide tissue's open conservation analysis for these region.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process
#' from other packages, such as esATAC.
#' @param bedInput \code{Character} scalar.
#' The directory of region BED file for analysis.
#' @param openConserveBedInput \code{Character} scalar.
#' The open level BED file for analysis. The first three columns are chromosome, start and end,
#' The remaining columns are region name and conservation score.
#' @param bedOutput \code{Character} scalar.
#' The BED output file directory of merged BED files.
#' Default: NULL (generated base on bedInput)
#' @param distrPdfOutput \code{Character} scalar.
#' The open conservation distribution figure for each tissue will be provided in PDF file.
#' @param ... Additional arguments, currently unused.
#' @details
#' We collected 201 DNase-seq or ATAC-seq sample from ENCODE and calculate their open level value.
#' They can be download and install automatically. So users do not need to configure themselves.
#' @return An invisible \code{\link{EnrichStep-class}}
#' object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{unzipAndMergeBed}}


#' @examples
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed.gz")
#' rs <- unzipAndMergeBed(bedInput = foregroundBedPath) %>%
#'          enrichTissueOpennessSpecificity


setGeneric("enrichTissueOpennessConserve",
           function(prevStep,
                    bedInput = NULL,
                    openConserveBedInput = NULL,
                    bedOutput = NULL,
                    distrPdfOutput = NULL,
                    ...) standardGeneric("enrichTissueOpennessConserve"))



#' @rdname TissueOpennessSpecificity
#' @aliases enrichTissueOpennessConserve
#' @export
setMethod(
    f = "enrichTissueOpennessConserve",
    signature = "Step",
    definition = function(prevStep,
                          bedInput = NULL,
                          openConserveBedInput = NULL,
                          bedOutput = NULL,
                          distrPdfOutput = NULL, ...){
        allpara <- c(list(Class = "TissueOpennessConserve",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname TissueOpennessConserve
#' @aliases tissueOpennessConserve
#' @export
tissueOpennessConserve <- function(bedInput,
                                   openConserveBedInput = NULL,
                                   bedOutput = NULL,
                                   distrPdfOutput = NULL, ...){
    allpara <- c(list(Class = "TissueOpennessConserve",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
