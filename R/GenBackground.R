setClass(Class = "GenBackground",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "GenBackground",
    definition = function(.Object,prevSteps,...){
        allparam <- list(...)
        message("aaaaaaaaa")
        print(allparam)
        inputForgroundBed <- allparam[["inputForgroundBed"]]
        bsgenome <- allparam[["bsgenome"]]
        outputForgroundBed <- allparam[["outputForgroundBed"]]
        outputBackgroundBed <- allparam[["outputBackgroundBed"]]
        outputAllRegionBed <- allparam[["outputAllRegionBed"]]
        rangeLen <- allparam[["rangeLen"]]
        sampleNumb <- allparam[["sampleNumb"]]
        message("ttttttttttttt")
        print(prevSteps)
        message("ttttttttttttt")
        if(!is.null(prevSteps)){
            forgroundBed <- getParam(prevSteps,"bedOutput")
            .Object@inputList[["inputForgroundBed"]] <- forgroundBed
        }

        if(!is.null(inputForgroundBed)){
            .Object@inputList[["inputForgroundBed"]] <- inputForgroundBed
        }

        if(is.null(outputForgroundBed)){
            .Object@outputList[["outputForgroundBed"]] <- getAutoPath(.Object,originPath = inputForgroundBed,regexProcName = "bed",suffix = "forground.bed")
        }else{
            .Object@outputList[["outputForgroundBed"]] <- outputForgroundBed
        }

        if(is.null(outputBackgroundBed)){
            .Object@outputList[["outputBackgroundBed"]] <- getAutoPath(.Object,originPath = inputForgroundBed,regexProcName = "bed",suffix = "background.bed")
        }else{
            .Object@outputList[["outputBackgroundBed"]] <- outputBackgroundBed
        }

        if(is.null(outputAllRegionBed)){
            .Object@outputList[["outputAllRegionBed"]] <- getAutoPath(.Object,originPath = inputForgroundBed,regexProcName = "bed",suffix = "allregion.bed")
        }else{
            .Object@outputList[["outputAllRegionBed"]] <- outputAllRegionBed
        }
        if(is.null(bsgenome)){
            .Object@paramList[["bsgenome"]] <- getRefRc("bsgenome")
        }else{
            stopifnot(is(bsgenom,"BSgenome"))
            .Object@paramList[["bsgenome"]] <- getRefRc("bsgenome")
        }
        .Object@paramList[["rangeLen"]] <- rangeLen
        message("rangeLen")
        print(.Object@paramList[["rangeLen"]] )
        .Object@paramList[["sampleNumb"]] <- sampleNumb
        .Object
    }
)


randomSampleOnGenome<-function(rangeLen, sampleNumber,bsgenome){
    chrlens <-seqlengths(bsgenome)
    selchr <- grep("_|M",names(chrlens),invert=TRUE)
    chrlens <- chrlens[selchr]
    startchrlens <- chrlens - rangeLen
    totallen <- sum(startchrlens)
    spnb <- floor(runif(sampleNumber) * totallen) + 1
    acclen <- startchrlens
    for(i in 2:length(acclen)){
        acclen[i] <- acclen[i-1] + acclen[i]
    }
    acclen <- c(0,acclen)
    gr <- GRanges()
    for(i in 1:(length(acclen)-1)){
        sel <- spnb[(spnb>acclen[i]) & (spnb<=acclen[i+1])]
        sel <- sel - acclen[i]
        gr <- c(gr,GRanges(seqnames = names(acclen)[i+1], ranges = IRanges(start = sel, width = 1000)))
    }
    return(sort(gr,ignore.strand=TRUE))
}

setMethod(
    f = "processing",
    signature = "GenBackground",
    definition = function(.Object,...){
        inputForgroundBed <- getParam(.Object,"inputForgroundBed")
        bsgenome <- getParam(.Object,"bsgenome")
        outputForgroundBed <- getParam(.Object,"outputForgroundBed")
        outputBackgroundBed <- getParam(.Object,"outputBackgroundBed")
        outputAllRegionBed <- getParam(.Object,"outputAllRegionBed")
        rangeLen <- getParam(.Object,"rangeLen")
        sampleNumb <- getParam(.Object,"sampleNumb")

        print(.Object@inputList)
        print(.Object@outputList)
        print(.Object@paramList)


        message("aaaaa")
        print(inputForgroundBed)
        print(rangeLen)
        foregroundgr <- import(con = inputForgroundBed,format = "bed")
        midpoint <- (start(foregroundgr) + end(foregroundgr))/2
        start(foregroundgr) <- floor(midpoint - rangeLen/2)
        end(foregroundgr) <- floor(midpoint + rangeLen/2)
        foregroundgr <- sort(foregroundgr,ignore.strand=TRUE)
        export.bed(object = foregroundgr, con = outputForgroundBed)
        backgroundgr <- randomSampleOnGenome(rangeLen, sampleNumb, bsgenome)
        export.bed(object = backgroundgr, con = outputBackgroundBed)
        regiongr <- c(foregroundgr,backgroundgr)
        export.bed(object = regiongr, con = outputAllRegionBed)
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "GenBackground",
    definition = function(.Object,...){
        if(is.null(.Object@inputList[["inputForgroundBed"]])){
            stop("inputForgroundBed is required.")
        }

    }
)



setMethod(
    f = "checkAllPath",
    signature = "GenBackground",
    definition = function(.Object,...){
        checkFileExist(.Object@inputList[["inputForgroundBed"]]);

    }
)

setMethod(
    f = "getReportValImp",
    signature = "GenBackground",
    definition = function(.Object,item,...){
        txt <- readLines(.Object@paramlist[["reportOutput"]])
        if(item == "total"){
            s<-strsplit(txt[1]," ")
            return(as.integer(s[[1]][1]))
        }
        if(item == "maprate"){
            s<-strsplit(txt[length(txt)],"% ")
            return(as.numeric(s[[1]][1])/100)
        }
        if(item == "detail"){
            return(txt)
        }
        stop(paste0(item," is not an item of report value."))
    }
)

setMethod(
    f = "getReportItemsImp",
    signature = "GenBackground",
    definition = function(.Object, ...){
        return(c("total","maprate","detail"))
    }
)

#' @name Bowtie2Mapping
#' @importFrom Rbowtie2 bowtie2
#' @title Use bowtie2 aligner to map reads to reference genome
#' @description
#' Use bowtie2 aligner to map reads to reference genome
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacRemoveAdapter}}
#' \code{\link{removeAdapter}}
#' @param reportOutput \code{Character} scalar.
#' The prefix of report files path.
#' @param bt2Idx \code{Character} scalar.
#' bowtie2 index files
#' prefix: 'dir/basename'
#' (minus trailing '.*.bt2' of 'dir/basename.*.bt2').
#' @param samOutput \code{Character} scalar.
#' A path to a SAM file
#' used for the alignment output.
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates
#' paired with file paths in fastqInput2.
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param fastqInput2 \code{Character} vector. It contains file paths with
#' #2 mates paired with file paths in fastqInput1.
#' For single-end sequencing files and interleaved paired-end
#' sequencing files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param paramList Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param threads \code{Integer} scalar.
#' The threads will be created in this process. default: 1
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bowtie2Mapping} instead.
#’ All additional arguments in paramList are interpreted as
#' additional parameters to be passed on to
#' bowtie2. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --no-mixed")
#' with white space splited just like command line,
#' or put them as \code{Character} vector
#' (e.g. c("--threads","8","--no-mixed")). Note that some
#' arguments("-x","--interleaved","-U","-1","-2","-S","threads") to the
#' bowtie2 are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{bowtie2_usage()} for details about available parameters.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacRemoveAdapter}}
#' \code{\link{removeAdapter}}
#' \code{\link{bowtie2}}
#' \code{\link{bowtie2_build}}
#' \code{\link{bowtie2_usage}}
#' \code{\link{atacSam2Bam}}
#' \code{\link{atacSamToBed}}
#' \code{\link{atacLibComplexQC}}
#' @examples
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' ## Building a bowtie2 index
#' library("Rbowtie2")
#' refs <- dir(system.file(package="esATAC", "extdata", "bt2","refs"),
#' full=TRUE)
#' bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' ## Alignments
#' reads_1 <- system.file(package="esATAC", "extdata", "bt2", "reads",
#' "reads_1.fastq")
#' reads_2 <- system.file(package="esATAC", "extdata", "bt2", "reads",
#' "reads_2.fastq")
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     (bowtie2Mapping(NULL,bt2Idx = file.path(td, "lambda_virus"),
#'        samOutput = file.path(td, "result.sam"),
#'        fastqInput1=reads_1,fastqInput2=reads_2,threads=3))
#'     head(readLines(file.path(td, "result.sam")))
#' }



setGeneric("enrichGenBackground",function(prevStep,
                                          inputForgroundBed = NULL,
                                          genome = NULL,
                                          outputForgroundBed = NULL,
                                          outputBackgroundBed = NULL,
                                          outputAllRegionBed = NULL,
                                          rangeLen = 1000,
                                          sampleNumb = 10000,
                                          ...) standardGeneric("enrichGenBackground"))



#' @rdname Bowtie2Mapping
#' @aliases atacBowtie2Mapping
#' @export
setMethod(
    f = "enrichGenBackground",
    signature = "Step",
    definition = function(prevStep,
                          inputForgroundBed = NULL,
                          genome = NULL,
                          outputForgroundBed = NULL,
                          outputBackgroundBed = NULL,
                          outputAllRegionBed = NULL,
                          rangeLen = 1000,
                          sampleNumb = 10000,
                          ...){
        allpara <- c(as.list(environment()),list(...))
        allpara[["Class"]] <- "GenBackground"
        allpara[["prevSteps"]] <- NULL
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname GenBackground
#' @aliases genBackground
#' @export
genBackground <- function(inputForgroundBed,
                          genome = NULL,
                          outputForgroundBed = NULL,
                          outputBackgroundBed = NULL,
                          outputAllRegionBed = NULL,
                          rangeLen = 1000,
                          sampleNumb = 10000,
                          ...){
    allpara <- c(list(Class = "GenBackground", prevSteps = NULL),as.list(environment()),list(...))
    print(allpara)
    step <- do.call(new,allpara)
    invisible(step)
}
