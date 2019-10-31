
setClass(Class = "SingleSampleReport",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "SingleSampleReport",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        htmlOutput <- allparam[["htmlOutput"]]


        if(is.null(htmlOutput)){
            output(.Object)$htmlOutput <- getStepWorkDir(.Object, filename = "report.html")
        }else{
            output(.Object)$htmlOutput <- htmlOutput
        }

        .Object
    }
)

#' @importFrom rmarkdown render

setMethod(
    f = "processing",
    signature = "SingleSampleReport",
    definition = function(.Object, ...){
        htmlOutput <- getParam(.Object, "htmlOutput")
        prevSteps <- list(...)[["prevSteps"]]
        prevStepsType <- lapply(prevSteps, function(step){
            return(stepType(step))
        })
        names(prevSteps) <- unlist(prevStepsType)

        save(prevSteps, file = getStepWorkDir(.Object = .Object, filename = "PrevSteps.Rdata"))

        reportmkd <- getStepWorkDir(.Object = .Object, filename = "Report.Rmd")

        file.copy(from = system.file(package = "enrichTF", "extdata","Report.Rmd"),
                  to = reportmkd)

        render(reportmkd)


        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "SingleSampleReport",
    definition = function(.Object,...){

    }
)



setMethod(
    f = "checkAllPath",
    signature = "SingleSampleReport",
    definition = function(.Object,...){


    }
)


#' @name SingleSampleReport
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Tissue's open conservation of the given region
#' @description
#' User provide region through a BED file.
#' This function will provide tissue's open conservation analysis for these region.
#' @param htmlOutput \code{Character} scalar.
#' The BED output file directory of merged BED files.
#' Default: NULL (generated base on bedInput)
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


setGeneric("enrichSingleSampleReport",
           function(prevStep, htmlOutput = NULL,...)
               standardGeneric("enrichSingleSampleReport"))



#' @rdname SingleSampleReport
#' @aliases enrichSingleSampleReport
#' @export
setMethod(
    f = "enrichSingleSampleReport",
    signature = "Step",
    definition = function(prevStep, htmlOutput = NULL, ...){
        allpara <- c(list(Class = "SingleSampleReport",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
