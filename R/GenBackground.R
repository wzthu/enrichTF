#' @importFrom methods new
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BiocGenerics start end
#' @importFrom BiocGenerics start<- end<-
setClass(Class = "GenBackground",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "GenBackground",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        inputForegroundBed <- allparam[["inputForegroundBed"]]
        genome <- allparam[["genome"]]
        outputForegroundBed <- allparam[["outputForegroundBed"]]
        outputBackgroundBed <- allparam[["outputBackgroundBed"]]
        outputRegionBed <- allparam[["outputRegionBed"]]
        regionLen <- allparam[["regionLen"]]
        sampleNumb <- allparam[["sampleNumb"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            foregroundBed <- getParam(prevStep,"bedOutput")
            input(.Object,"inputForegroundBed") <- foregroundBed
        }

        if(!is.null(inputForegroundBed)){
            input(.Object,"inputForegroundBed") <- inputForegroundBed
        }

        if(is.null(outputForegroundBed)){
            output(.Object,"outputForegroundBed") <- getAutoPath(.Object,originPath = input(.Object)[["inputForegroundBed"]],regexSuffixName = "bed",suffix = "foreground.bed")
        }else{
            output(.Object,"outputForegroundBed") <- outputForegroundBed
        }

        if(is.null(outputBackgroundBed)){
            output(.Object,"outputBackgroundBed") <- getAutoPath(.Object,originPath = input(.Object)[["inputForegroundBed"]],regexSuffixName = "bed",suffix = "background.bed")
        }else{
            output(.Object,"outputBackgroundBed") <- outputBackgroundBed
        }

        if(is.null(outputRegionBed)){
            output(.Object,"outputRegionBed") <- getAutoPath(.Object,originPath = input(.Object)[["inputForegroundBed"]],regexSuffixName = "bed",suffix = "allregion.bed")
        }else{
            output(.Object,"outputRegionBed") <- outputRegionBed
        }
        if(is.null(genome)){
            if(getGenome() == "testgenome"){
                param(.Object,"bsgenome") <- BSgenome::getBSgenome("hg19")
            }else{
                param(.Object,"bsgenome") <- BSgenome::getBSgenome(getGenome())
            }
        }else{
            if(genome == "testgenome"){
                genome <- "hg19"
            }
            param(.Object,"bsgenome") <- BSgenome::getBSgenome(genome = genome)


        }
        param(.Object,"regionLen") <- regionLen
        if(is.null(sampleNumb)){
            param(.Object,"sampleNumb")<- 0
        }else{
            param(.Object,"sampleNumb") <- sampleNumb
        }

        .Object
    }
)


# randomSampleOnGenome<-function(regionLen, sampleNumber,bsgenome){
#     chrlens <-seqlengths(bsgenome)
#     selchr <- grep("_|M",names(chrlens),invert=TRUE)
#     chrlens <- chrlens[selchr]
#     startchrlens <- chrlens - regionLen
#     totallen <- sum(startchrlens)
#     spnb <- floor(runif(sampleNumber) * totallen) + 1
#     acclen <- startchrlens
#     for(i in 2:length(acclen)){
#         acclen[i] <- acclen[i-1] + acclen[i]
#     }
#     acclen <- c(0,acclen)
#     gr <- GRanges()
#     for(i in 1:(length(acclen)-1)){
#         sel <- spnb[(spnb>acclen[i]) & (spnb<=acclen[i+1])]
#         sel <- sel - acclen[i]
#         gr <- c(gr,GRanges(seqnames = names(acclen)[i+1], ranges = IRanges(start = sel, width = 1000)))
#     }
#     return(sort(gr,ignore.strand=TRUE))
# }

randomSampleOnGenome<-function(regionLen, sampleNumber,bsgenome){
    chrlens <-seqlengths(bsgenome)
    if(getGenome() == "testgenome"){
        chrlens<-chrlens['chr1']
    }
    selchr <- grep("_|M",names(chrlens),invert=TRUE)
    chrlens <- chrlens[selchr]
    startchrlens <- chrlens - regionLen
    spchrs <- sample(x = names(startchrlens),size =  sampleNumber, replace = TRUE, prob = startchrlens / sum(startchrlens))
    gr <- GRanges()
    for(chr in names(startchrlens)){
        if(sum(spchrs == chr) == 0){
            next
        }
        startpt <- sample(x = 1:startchrlens[chr],size = sum(spchrs == chr),replace = FALSE)
        #non overlapped method:
        #startpt <- sample(x = 1:(startchrlens[chr] - sum(spchrs == chr) * regionLen),size = sum(spchrs == chr),replace = FALSE)
        #startpt <- startpt + 0:(sum(spchrs == chr)-1) * regionLen
        gr <- c(gr,GRanges(seqnames = chr, ranges = IRanges(start = startpt, width = 1000)))
    }
    return(sort(gr,ignore.strand=TRUE))
}

setMethod(
    f = "processing",
    signature = "GenBackground",
    definition = function(.Object,...){
        inputForegroundBed <- getParam(.Object,"inputForegroundBed")
        bsgenome <- getParam(.Object,"bsgenome")
        outputForegroundBed <- getParam(.Object,"outputForegroundBed")
        outputBackgroundBed <- getParam(.Object,"outputBackgroundBed")
        outputRegionBed <- getParam(.Object,"outputRegionBed")
        regionLen <- getParam(.Object,"regionLen")
        sampleNumb <- getParam(.Object,"sampleNumb")

        foregroundgr <- import(con = inputForegroundBed,format = "bed")
        midpoint <- (start(foregroundgr) + end(foregroundgr))/2
        start(foregroundgr) <- floor(midpoint - regionLen/2)
        end(foregroundgr) <- floor(midpoint + regionLen/2)
        foregroundgr <- sort(foregroundgr,ignore.strand=TRUE)
        mcols(foregroundgr)$name <-  1:length(foregroundgr)
        export.bed(object = foregroundgr, con = outputForegroundBed)
        if(sampleNumb == 0){
            sampleNumb = length(foregroundgr)
        }
        backgroundgr <- randomSampleOnGenome(regionLen, sampleNumb, bsgenome)
        mcols(backgroundgr)$name <- (length(foregroundgr) + 1) : (length(foregroundgr) + length(backgroundgr))
        export.bed(object = backgroundgr, con = outputBackgroundBed)
        regiongr <- c(foregroundgr,backgroundgr)
        export.bed(object = regiongr, con = outputRegionBed)
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "GenBackground",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["inputForegroundBed"]])){
            stop("inputForegroundBed is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "GenBackground",
    definition = function(.Object,...){
        checkFileExist(input(.Object)[["inputForegroundBed"]]);

    }
)

setMethod(
    f = "getReportValImp",
    signature = "GenBackground",
    definition = function(.Object,item,...){
        txt <- readLines(param(.Object)[["reportOutput"]])
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


#' @name GenBackground
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Generate background regions and reset the size of foreground regions
#' @description
#' Use uniform distribution to generate background sequence regions from genome.
#' The size of foreground regions will be unified into the length specified in argument.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process from other packages, such as esATAC.
#' @param inputForegroundBed \code{Character} scalar.
#' The directory of foreground BED file.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome such as "hg19", "mm10", etc.
#' Default: NULL (e.g. after library (enrichTF), you can call function \code{setGenome("hg19")})
#' @param outputForegroundBed \code{Character} scalar.
#' The BED file directory of reshaped foreground regions.
#' Default: NULL (generated base on inputForegroundBed)
#' @param outputBackgroundBed \code{Character} scalar.
#' The BED file directory of reshaped background regions.
#' Default: NULL (generated base on inputForegroundBed)
#' @param outputRegionBed \code{Character} scalar.
#' Foreground and background merged BED files.
#' Default: NULL (generated base on inputForegroundBed)
#' @param regionLen \code{Character} scalar. It sets the length of forground sequence regions. Default: 1000
#' @param sampleNumb \code{numeric} scalar.
#' It sets the number of background regions that will be sampled. Default: 10000
#' @param ... Additional arguments, currently unused.
#' @details
#' Use uniform distribution to generate background sequence regions from genome.
#' The size of foreground regions will be unified into the length specified in argument.
#' @return An invisible \code{\link{EnrichStep-class}} object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{regionConnectTargetGene}}
#' \code{\link{findMotifsInRegions}}
#' \code{\link{tfsEnrichInRegions}}

#' @examples
#' setGenome("testgenome") #Use "hg19","hg38",etc. for your application
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' gen <- genBackground(inputForegroundBed = foregroundBedPath)


setGeneric("enrichGenBackground",function(prevStep,
                                          inputForegroundBed = NULL,
                                          genome = NULL,
                                          outputForegroundBed = NULL,
                                          outputBackgroundBed = NULL,
                                          outputRegionBed = NULL,
                                          regionLen = 1000,
                                          sampleNumb = 10000,
                                          ...) standardGeneric("enrichGenBackground"))



#' @rdname GenBackground
#' @aliases enrichGenBackground
#' @export
setMethod(
    f = "enrichGenBackground",
    signature = "Step",
    definition = function(prevStep,
                          inputForegroundBed = NULL,
                          genome = NULL,
                          outputForegroundBed = NULL,
                          outputBackgroundBed = NULL,
                          outputRegionBed = NULL,
                          regionLen = 1000,
                          sampleNumb = NULL,
                          ...){
        allpara <- c(list(Class = "GenBackground", prevSteps = list(prevStep)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname GenBackground
#' @aliases genBackground
#' @export
genBackground <- function(inputForegroundBed,
                          genome = NULL,
                          outputForegroundBed = NULL,
                          outputBackgroundBed = NULL,
                          outputRegionBed = NULL,
                          regionLen = 1000,
                          sampleNumb = NULL,
                          ...){
    allpara <- c(list(Class = "GenBackground", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
