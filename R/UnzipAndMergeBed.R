
setClass(Class = "UnzipAndMergeBed",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "UnzipAndMergeBed",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        bedOutput <- allparam[["bedOutput"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            bedInput0 <- getParam(prevStep,"bedOutput")
            input(.Object)$bedInput <- bedInput0
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



#' @importFrom R.utils isBzipped
#' @importFrom R.utils isGzipped
#' @importFrom R.utils bunzip2
#' @importFrom R.utils gunzip

setGeneric(
    name = "decompressBed",
    def = function(.Object,filename,destpath,...){
        standardGeneric("decompressBed")
    }
)
setMethod(
    f = "decompressBed",
    signature = "UnzipAndMergeBed",
    definition = function(.Object,filename,destpath,...){
        destname<-file.path(destpath,basename(filename))
        writeLog(.Object,paste0("processing file:"))
        writeLog(.Object,sprintf("source:%s",filename))
        writeLog(.Object,sprintf("destination:%s",destname))
        if(isBzipped(filename)){
            destname<-gsub(sprintf("[.]%s$", "bz2"), "", destname, ignore.case=TRUE)
            return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else if(isGzipped(filename)){
            destname<-gsub(sprintf("[.]%s$", "gz"), "", destname, ignore.case=TRUE)
            return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else{
            return(filename)
        }
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
            tmpfile <- decompressBed(.Object, bedfile, getStepWorkDir(.Object))
            bed <- read.table(tmpfile, header = FALSE, sep = "\t")
            bed <- bed[,1:3]
            colnames(bed) <-c("chrom", "start", "end")
            if(bedfile != tmpfile){
                unlink(tmpfile)
            }
            return(bed)
        })

        bed <- do.call(rbind, args = beds)
        write.table(bed, file = bedOutput, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

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
#' @title Unzip all zipped BED files and
#' merge them into one BED file
#' @description
#' This function process region BED files in three step:
#' First, unzip the gzip and bzip2  BED input files.
#' Second, select first 3 columns of the BED files.
#' Third, merge the BED files into one BED.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process
#' from other packages, such as ATAC-seq
#' peak calling result from esATAC.
#' @param bedInput \code{Character} scalar or vector.
#' The directory of region BED files for analysis.
#' BED, BED.gz, BED.bz2 formats are supported.
#' @param bedOutput \code{Character} scalar.
#' The BED output file directory of merged BED files.
#' Default: NULL (generated base on first BED file in  bedInput)
#' @param ... Additional arguments, currently unused.
#' @details
#' All compressed files will be de-compressed.
#' Only first 3 columns (chromasomes, start and end) will be collected.
#' All BED files will be merged into one BED file.
#' @return An invisible \code{\link{EnrichStep-class}}
#' object (inherit from \code{\link{Step-class}}) scalar for downstream analysis.
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



#' @rdname UnzipAndMergeBed
#' @aliases enrichUnzipAndMergeBed
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
#' @rdname UnzipAndMergeBed
#' @aliases unzipAndMergeBed
#' @export
unzipAndMergeBed <- function(bedInput, bedOutput = NULL, ...){
    allpara <- c(list(Class = "UnzipAndMergeBed",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
