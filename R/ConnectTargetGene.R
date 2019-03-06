setClass(Class = "RegionConnectTargetGene",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        print(allparam)
        inputForegroundBed <- allparam[["inputForegroundBed"]]
        inputBackgroundBed <- allparam[["inputBackgroundBed"]]
        outputForegroundBed <- allparam[["outputForegroundBed"]]
        outputBackgroundBed <- allparam[["outputBackgroundBed"]]
        regularGeneCorrBed <- allparam[["regularGeneCorrBed"]]
        enhancerRegularGeneCorrBed <- allparam[["enhancerRegularGeneCorrBed"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            .Object@inputList[["inputForegroundBed"]] <- getParam(prevStep,"outputForegroundBed")
            .Object@inputList[["inputBackgroundBed"]] <- getParam(prevStep,"outputBackgroundBed")
        }


        if(is.null(outputForegroundBed)){
            .Object@outputList[["outputForegroundBed"]] <- getAutoPath(.Object,originPath = .Object@inputList[["inputForegroundBed"]],regexProcName = "foreground.bed",suffix = "gene.foreground.bed")
        }else{
            .Object@outputList[["outputForegroundBed"]] <- outputForegroundBed
        }


        if(is.null(outputBackgroundBed)){
            .Object@outputList[["outputBackgroundBed"]] <- getAutoPath(.Object,originPath = .Object@inputList[["inputBackgroundBed"]],regexProcName = "background.bed",suffix = "gene.background.bed")
        }else{
            .Object@outputList[["outputBackgroundBed"]] <- outputBackgroundBed
        }

        if(is.null(regularGeneCorrBed)){
            .Object@outputList[["regularGeneCorrBed"]] <- getRefFiles("RE_gene_corr")
        }else{
            .Object@outputList[["regularGeneCorrBed"]] <- regularGeneCorrBed
        }

        if(is.null(enhancerRegularGeneCorrBed)){
            .Object@outputList[["enhancerRegularGeneCorrBed"]] <- getRefFiles("Enhancer_RE_gene_corr")
        }else{
            .Object@outputList[["enhancerRegularGeneCorrBed"]] <- enhancerRegularGeneCorrBed
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
        enhancerRegularGeneCorrBed <- getParam(.Object,"enhancerRegularGeneCorrBed")

        print(.Object@inputList)
        print(.Object@outputList)
        print(.Object@paramList)

        inputForegroundgr <- import(con=inputForegroundBed)
        inputBackgroundgr <- import(con=inputBackgroundBed)

        rg <- import(con=regularGeneCorrBed,colnames=c("name","score","blockCount"))
        erg <- import(con=regularGeneCorrBed,colnames=c("name","score","blockCount"))

        pairs <- findOverlapPairs(inputForegroundgr,rg)
        first(pairs)$name  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputForegroundgr <- first(pairs)

        pairs <- findOverlapPairs(inputForegroundgr,erg)
        first(pairs)$name  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputForegroundgr <- c(outputForegroundgr,first(pairs))

        pairs <- findOverlapPairs(inputBackgroundgr,rg)
        first(pairs)$name  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputBackgroundgr <- first(pairs)

        pairs <- findOverlapPairs(inputBackgroundgr,erg)
        first(pairs)$name  <- second(pairs)$name
        first(pairs)$score  <- second(pairs)$score
        first(pairs)$blockCount  <- second(pairs)$blockCount

        outputBackgroundgr <- c(outputBackgroundgr,first(pairs))

        write.table(as.data.frame(outputForegroundgr)[,c("seqnames","start","end","name","score",
                                                         "strand","strand","strand","strand","blockCount")],
                    outputForegroundBed,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
        write.table(as.data.frame(outputBackgroundgr)[,c("seqnames","start","end","name","score",
                                                         "strand","strand","strand","strand","blockCount")],
                    outputBackgroundBed,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
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
        if(is.null(.Object@inputList[["inputForegroundBed"]])){
            stop("inputForegroundBed is required.")
        }
        if(is.null(.Object@inputList[["inputBackgroundBed"]])){
            stop("inputBackgroundBed is required.")
        }

    }
)



setMethod(
    f = "checkAllPath",
    signature = "RegionConnectTargetGene",
    definition = function(.Object,...){
        checkFileExist(.Object@inputList[["inputForegroundBed"]])
        checkFileExist(.Object@inputList[["inputBackgroundBed"]])

    }
)

setMethod(
    f = "getReportValImp",
    signature = "RegionConnectTargetGene",
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
    signature = "RegionConnectTargetGene",
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



setGeneric("enrichRegionConnectTargetGene",function(prevStep,
                                                    inputForegroundBed = NULL,
                                                    inputBackgroundBed = NULL,
                                                    outputForegroundBed = NULL,
                                                    outputBackgroundBed = NULL,
                                                    regularGeneCorrBed = NULL,
                                                    enhancerRegularGeneCorrBed = NULL,
                                                    ...) standardGeneric("enrichRegionConnectTargetGene"))



#' @rdname MotifsInRegions
#' @aliases enrichMotifsInRegions
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
                          ...){
        allpara <- c(list(Class = "RegionConnectTargetGene", prevSteps = list(prevStep)),as.list(environment()),list(...))
        print(allpara)
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname MotifsInRegions
#' @aliases motifsInRegions
#' @export
regionConnectTargetGene <- function(inputForegroundBed,
                                    inputBackgroundBed,
                                    outputForegroundBed = NULL,
                                    outputBackgroundBed = NULL,
                                    regularGeneCorrBed = NULL,
                                    enhancerRegularGeneCorrBed = NULL,
                                    ...){
    allpara <- c(list(Class = "RegionConnectTargetGene", prevSteps = list()),as.list(environment()),list(...))
    print(allpara)
    step <- do.call(new,allpara)
    invisible(step)
}
