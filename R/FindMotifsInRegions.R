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
            second(motif_region_pair)$motifName <-paste(pwmObj[[i]]@tags$motifName,pwmObj[[i]]@tags$motifSrc,pwmObj[[i]]@tags$motifPlf,sep = '/')
            if(i == 1){
                result <- second(motif_region_pair)
            }else{
                result <- c(result,second(motif_region_pair))
            }

        }

        .Object@propList[["motifs_in_region"]] <- result

        write.table(as.data.frame(result)[,c("seqnames","start","end","name","score","motifName")],
                    file = outputRegionMotifBed, sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

   #     export.bed(result,con = outputRegionMotifBed)

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



#' @name MotifsInRegions
#' @importFrom motifmatchr matchMotifs
#' @title Find motifs in all input sequence regions
#' @description
#' Scan for motif occurrences using the prepared PWMs and obtain the promising candidate motifs in these regions.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' It needs to be the return value of upstream process from \code{\link{genBackground}} or \code{\link{enrichGenBackground}}
#' @param inputRegionBed \code{Character} scalar.
#' BED file for regions including foreground and background sequences.
#' @param outputRegionMotifBed \code{Character} scalar.
#' BED file for regions with motif candidates.
#' Default: NULL (generated base on inputForegroundBed)
#' @param motifRc \code{Character} scalar.
#' Motif Resources can be one of "integrate" 
#' (integrated by us and can be download from internet automatically 
#' if call the function \code{setGenome("hg19")}), 
#' "jaspar" package JASPAR2018, 
#' or "pwmfile" (User defined PWM file. inputPwmFile is required).
#' @param inputPwmFile \code{Character} scalar.
#' when "pwmfile" is set for motifRc, use this argument to provide PWM file directory.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome, such as "hg19", "mm10", etc. 
#' Default: NULL (e.g. after \code{library (enrichTF)}, you can call function \code{setGenome("hg19")})
#' @param ... Additional arguments, currently unused.
#' @details
#' Scan for motif occurrences using the prepared PWMs and 
#' obtain the promising candidate motifs in these regions.
#' @return An invisible \code{\link{EnrichTF-class}} object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{genBackground}}
#' \code{\link{findMotifsInRegions}}
#' \code{\link{tfsEnrichInRegions}}
#' @examples
#' setGenome("hg19")
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testForeGround.bed")
#' gen <- genBackground(inputForegroundBed = foregroundBedPath)
#' findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate")





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
