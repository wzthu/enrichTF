setClass(Class = "FindMotifsInRegions",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "FindMotifsInRegions",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        print(allparam)
        inputRegionBed <- allparam[["inputRegionBed"]]
        motifRc <- allparam[["motifRc"]]
        inputPwmFile <- allparam[["inputPwmFile"]]
        genome <- allparam[["genome"]]
        outputRegionMotifBed <- allparam[["outputRegionMotifBed"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            outputRegionBed <- getParam(prevStep,"outputRegionBed")
            .Object@inputList[["inputRegionBed"]] <- outputRegionBed
        }
        motifRc <- match.arg(motifRc, c("integrate","jaspar","pwmfile"))

        .Object@paramList[["motifRc"]] <- motifRc

        if(is.null(outputRegionMotifBed)){
            .Object@outputList[["outputRegionMotifBed"]] <- getAutoPath(.Object,originPath = .Object@inputList[["inputRegionBed"]],regexProcName = "allregion.bed",suffix = "region.motif.bed")
        }else{
            .Object@outputList[["outputRegionMotifBed"]] <- outputRegionMotifBed
        }

        if(motifRc == "integrate"){
            load(getRefFiles("motifPWMOBJ"))
            .Object@inputList[["pwmObj"]] <- pwm
        }else if(motifRc == "jarspar"){
            .Object@inputList[["pwmObj"]] <- JASPAR2018::JASPAR2018

        }else if(motifRc == "pwmfile"){
            .Object@inputList[["pwmObj"]] <- inputPwmFile
        }



        if(is.null(genome)){
            .Object@paramList[["genome"]] <- getGenome()
        }else{
            .Object@paramList[["genome"]] <- genome
        }
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "FindMotifsInRegions",
    definition = function(.Object,...){
        inputRegionBed <- getParam(.Object,"inputRegionBed")
        pwmObj <- getParam(.Object,"pwmObj")
        genome <- getParam(.Object,"genome")
        outputRegionMotifBed <- getParam(.Object,"outputRegionMotifBed")
        regions <- import(con = inputRegionBed,format = "bed")
        print(.Object@inputList)
        print(.Object@outputList)
        print(.Object@paramList)
        motif_ix <-motifmatchr::matchMotifs(pwms = pwmObj, subject = regions, genome = genome, out="positions")
        result <- c()
        .Object@propList[["motif_ix"]] <-motif_ix
        for(i in 1:length(motif_ix)){
            motif_region_pair <- findOverlapPairs(motif_ix[[i]],regions,ignore.strand = TRUE)
            second(motif_region_pair)$score <- first(motif_region_pair)$score
            second(motif_region_pair)$name <-paste(pwmObj[[i]]@tags$motifName,pwmObj[[i]]@tags$motifSrc,pwmObj[[i]]@tags$motifPlf,sep = '/')
            if(i == 1){
                result <- second(motif_region_pair)
            }else{
                result <- c(result,second(motif_region_pair))
            }

        }

        .Object@propList[["motifs_in_region"]] <- result

        export.bed(result,con = outputRegionMotifBed)

        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "FindMotifsInRegions",
    definition = function(.Object,...){
        if(is.null(.Object@inputList[["inputRegionBed"]])){
            stop("inputRegionBed is required.")
        }

    }
)



setMethod(
    f = "checkAllPath",
    signature = "FindMotifsInRegions",
    definition = function(.Object,...){
        checkFileExist(.Object@inputList[["inputRegionBed"]])

    }
)

setMethod(
    f = "getReportValImp",
    signature = "FindMotifsInRegions",
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
    signature = "FindMotifsInRegions",
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



setGeneric("enrichFindMotifsInRegions",function(prevStep,
                                          inputRegionBed = NULL,
                                          outputRegionMotifBed = NULL,
                                          motifRc = c("integrate","jaspar","pwmfile"),
                                          inputPwmFile = NULL,
                                          genome = NULL,
                                          ...) standardGeneric("enrichFindMotifsInRegions"))



#' @rdname MotifsInRegions
#' @aliases enrichMotifsInRegions
#' @export
setMethod(
    f = "enrichFindMotifsInRegions",
    signature = "Step",
    definition = function(prevStep,
                          inputRegionBed = NULL,
                          outputRegionMotifBed = NULL,
                          motifRc = c("integrate","jaspar","pwmfile"),
                          inputPwmFile = NULL,
                          genome = NULL,
                          ...){
        allpara <- c(list(Class = "FindMotifsInRegions", prevSteps = list(prevStep)),as.list(environment()),list(...))
        print(allpara)
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname MotifsInRegions
#' @aliases motifsInRegions
#' @export
findMotifsInRegions <- function(inputRegionBed,
                                outputRegionMotifBed,
                                motifRc = c("integrate","jaspar","pwmfile"),
                                inputPwmFile = NULL,
                                genome = NULL,
                                ...){
    allpara <- c(list(Class = "FindMotifsInRegions", prevSteps = list()),as.list(environment()),list(...))
    print(allpara)
    step <- do.call(new,allpara)
    invisible(step)
}
