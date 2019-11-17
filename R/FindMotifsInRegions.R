#' @importFrom S4Vectors mcols<-
#' @importFrom S4Vectors mcols
#' @importFrom stats fisher.test p.adjust t.test
#' @importFrom S4Vectors second<- first<-
#' @importFrom rtracklayer export.bed
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand<-
#' @importFrom BiocGenerics strand
setClass(Class = "FindMotifsInRegions",
         contains = "EnrichStep"
)

setMethod(
    f = "init",
    signature = "FindMotifsInRegions",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        inputRegionBed <- allparam[["inputRegionBed"]]
        motifRc <- allparam[["motifRc"]]
        inputPwmFile <- allparam[["inputPwmFile"]]
        genome <- allparam[["genome"]]
        outputRegionMotifBed <- allparam[["outputRegionMotifBed"]]
        threads <- allparam[["threads"]]
        if(length(prevSteps)>0){
            prevStep <- prevSteps[[1]]
            outputRegionBed <- getParam(prevStep,"outputRegionBed")
            input(.Object)$inputRegionBed <- outputRegionBed
        }
        motifRc <- match.arg(motifRc, c("integrate","jaspar","pwmfile"))

        param(.Object)$motifRc <- motifRc

        if(!is.null(inputRegionBed)){
            input(.Object)$inputRegionBed <- inputRegionBed
        }

        if(is.null(outputRegionMotifBed)){
            output(.Object)$outputRegionMotifBed <-
                getAutoPath(.Object,originPath =
                                .Object$inputList[["inputRegionBed"]],
                            regexSuffixName = "allregion.bed",
                            suffix = "region.motif.bed")
        }else{
            output(.Object)$outputRegionMotifBed <- outputRegionMotifBed
        }

        output(.Object)$outputMotifBed <-
            getAutoPath(.Object,originPath =
                            .Object$inputList[["inputRegionBed"]],
                        regexSuffixName = "allregion.bed",
                        suffix = "motif.bed")

        if(motifRc == "integrate"){

            param(.Object)$pwmObj <- get(load(getRefFiles("motifPWMOBJ")))
        }else if(motifRc == "jarspar"){
            param(.Object)$pwmObj <- JASPAR2018::JASPAR2018

        }else if(motifRc == "pwmfile"){
            param(.Object)$pwmObj <- getMotifInfo1(inputPwmFile)
        }



        if(is.null(genome)){
            param(.Object)$genome <- getGenome()
        }else{
            param(.Object)$genome <- genome
        }
        if(.Object$paramList[["genome"]] == "testgenome"){
            param(.Object)$genome <- "hg19"
        }
        if(is.null(threads)){
            param(.Object)$threads <- getThread
        }else{
            param(.Object)$threads <- threads
        }
        .Object
    }
)


homerOutputToBed <- function(regionsbed, outputtxt){
    #  regionsbed <- "region.bed"
    #  outputtxt <- "outputfile.txt"

    regions <- import(con = regionsbed,format = "bed")
    txt <- read.table(outputtxt,sep = "\t", header = TRUE)
    txt[['idx']] <- 1:nrow(txt)
    start(regions) <- start(regions) -1

    pos <- txt[txt$Strand == '+',]
    neg <- txt[txt$Strand == '-',]

    mcols(regions)$mid <- round((start(regions)/2 + end(regions)/2))
    mcols(regions)$id <- 1:length(regions)

    posbed <- regions[pos$PositionID,]
    negbed <- regions[neg$PositionID,]

    posbedst <- mcols(posbed)$mid + as.integer(pos$Offset)

    posbed <- data.frame(chrom = seqnames(posbed),
                         start = posbedst,
                         end = posbedst + as.numeric(sapply(as.character(pos$Sequence),nchar)),
                         strand = '+',
                         names = as.character(pos$Motif.Name),
                         score = as.numeric(pos$MotifScore),
                         idx = as.numeric(pos$idx)
                         #           str = as.character(pos$Sequence)
    )



    negbeded <- mcols(negbed)$mid + as.integer(neg$Offset)

    negbed <- data.frame(chrom = seqnames(negbed),
                         start = negbeded - as.numeric(sapply(as.character(neg$Sequence),nchar)),
                         end = negbeded,
                         strand = '-',
                         names = as.character(neg$Motif.Name),
                         score = as.numeric(neg$MotifScore),
                         idx = as.numeric(neg$idx)
                         #                     str = as.character(neg$Sequence)
    )



    allbed <- rbind(posbed,negbed)

    idx <- order(allbed$idx)
    allbed <- allbed[idx,]

    allbed <- GRanges(allbed)
    names(allbed) <- mcols(allbed)$names
    return(allbed)
    #  export.bed(con="tobed.r.test.bed",allbed)

}


#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel makeCluster



setMethod(
    f = "processing",
    signature = "FindMotifsInRegions",
    definition = function(.Object,...){
        #        return(.Object)
        inputRegionBed <- getParam(.Object,"inputRegionBed")
        pwmObj <- getParam(.Object,"pwmObj")
        genome <- getParam(.Object,"genome")
        threads <- getParam(.Object,"threads")
        outputRegionMotifBed <- getParam(.Object,"outputRegionMotifBed")
        outputMotifBed <- getParam(.Object,"outputMotifBed")
        if(genome == "testgenome"){
            pwmObj = pwmObj[names(pwmObj)[(seq_len(length(pwmObj))%%4==0)]]
        }
        regions <- import(con = inputRegionBed,format = "bed")
        homer <- getRefFiles("HOMER")
        if(dir.exists(file.path(homer,"bin"))){
            homerinput <- file.path(getStepWorkDir(.Object), "homer.input.bed")
            system(paste("awk -F'\t' '{print $1\"\t\"$2\"\t\"$3}'",inputRegionBed, ">",homerinput))
            findMotifsGenome <- file.path(homer,"bin","findMotifsGenome.pl")
            homeroutput <- file.path(getStepWorkDir(.Object), "homer.output.txt")
            stopifnot(0==system(paste(findMotifsGenome, homerinput, getGenome(),
                                      file.path(getStepWorkDir(.Object),"region.all.hommer/"),
                                      "-find", getRefFiles("motifpwm"),
                                      ">",homeroutput)))
            motifranges <- homerOutputToBed(inputRegionBed,homeroutput)

            motif_region_pair <- findOverlapPairs(motifranges,
                                                  regions,ignore.strand = FALSE)
            if(length(second(motif_region_pair))>0){
                #            second(motif_region_pair)$score <-
                #                first(motif_region_pair)$score
                #            strand(second(motif_region_pair))<-
                #                strand(first(motif_region_pair))
                #            second(motif_region_pair)$motifName <- first(motif_region_pair)$name
                first(motif_region_pair)$name <- second(motif_region_pair)$name
                result <- second(motif_region_pair)
                result$score <- first(motif_region_pair)$score
                strand(result) <-strand(first(motif_region_pair))
                mcols(result)$motifName <- first(motif_region_pair)$names


                write.table(as.data.frame(result)[,c("seqnames","start","end",
                                                     "name","score","strand","motifName")],
                            file = outputRegionMotifBed, sep="\t",quote = FALSE,
                            row.names = FALSE,col.names = FALSE)

                result <- first(motif_region_pair)
                names(result) <- 1:length(result)

                write.table(as.data.frame(result)[,c("seqnames","start","end",
                                                     "names","score","strand","name")],
                            file = outputMotifBed, sep="\t",quote = FALSE,
                            row.names = FALSE,col.names = FALSE)
            }

        }else{
            motif_ix <- NULL
            if(osname == "osx" || osname == "linux"){
                motif_ix <-
                    parallel::mclapply(pwmObj,
                                       motifmatchr::matchMotifs,
                                       subject = regions,
                                       genome = genome,
                                       out = "positions", p.cutoff = 5e-04, mc.cores = threads)
            }else{
                cl <- makeCluster(threads)
                motif_ix <- parallel::parLapply(cl = cl, X = pwmObj,
                                   fun = motifmatchr::matchMotifs,
                                   subject = regions,
                                   genome = genome,
                                   out = "positions", p.cutoff = 5e-04)
                stopCluster(cl)

            }
            #motifmatchr::matchMotifs(pwms = pwmObj, subject = regions, genome = genome, out="positions")
            result <- c()
            #        .Object@propList[["motif_ix"]] <-motif_ix
            # for(i in 1:length(motif_ix)){
            #     motif_region_pair <- findOverlapPairs(motif_ix[[i]][[1]],regions,ignore.strand = TRUE)
            #     if(length(second(motif_region_pair))>0){
            #         second(motif_region_pair)$score <- first(motif_region_pair)$score
            #         second(motif_region_pair)$motifName <-pwmObj[[i]]@name
            #     }else{
            #         next
            #     }
            #     if(i == 1){
            #         result <- second(motif_region_pair)[second(motif_region_pair)$score >= pwmObj[[i]]@tags$threshold]
            #     }else{
            #         result <- c(result,second(motif_region_pair)[second(motif_region_pair)$score >= pwmObj[[i]]@tags$threshold])
            #     }
            #
            # }
            #print(head(motif_ix[[1]][[1]]))
            #print(findOverlapPairs(motif_ix[[1]][[1]],
            #                       regions,ignore.strand = FALSE))

            lapply(seq_len(length(motif_ix)), function(i){
                motif_region_pair <- findOverlapPairs(motif_ix[[i]][[1]],
                                                      regions,ignore.strand = FALSE)
                if(length(second(motif_region_pair))>0){
                    second(motif_region_pair)$score <-
                        first(motif_region_pair)$score
                    strand(second(motif_region_pair))<-
                        strand(first(motif_region_pair))
                    second(motif_region_pair)$motifName <- pwmObj[[i]]@name
                    first(motif_region_pair)$motifName <- pwmObj[[i]]@name
                    first(motif_region_pair)$name <- second(motif_region_pair)$name
                    result <- second(motif_region_pair)[
                        second(motif_region_pair)$score >=
                            pwmObj[[i]]@tags$threshold]

                    write.table(as.data.frame(result)[,c("seqnames","start","end",
                                                         "name","score","strand","motifName")],
                                file = outputRegionMotifBed, sep="\t",quote = FALSE, append = TRUE,
                                row.names = FALSE,col.names = FALSE)

                    result <- first(motif_region_pair)[
                        first(motif_region_pair)$score >=
                            pwmObj[[i]]@tags$threshold]

                    write.table(as.data.frame(result)[,c("seqnames","start","end",
                                                         "motifName","score","strand","name")],
                                file = outputMotifBed, sep="\t",quote = FALSE, append = TRUE,
                                row.names = FALSE,col.names = FALSE)
                }
                return(NULL)

                # return(second(motif_region_pair)[
                #        second(motif_region_pair)$score >=
                #           pwmObj[[i]]@tags$threshold])
            })
        }
        #  result <- do.call("c",result)
        #        .Object@propList[["motifs_in_region"]] <- result



        #     export.bed(result,con = outputRegionMotifBed)

        .Object
    }
)


#' @name MotifsInRegions
#' @importFrom motifmatchr matchMotifs
#' @importFrom JASPAR2018 JASPAR2018
#' @title Find motifs in all input sequence regions
#' @description
#' Scan for motif occurrences using the prepared PWMs and obtain
#' the promising candidate motifs in these regions.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' This parameter is available when the upstream step function
#' (printMap() to see the previous functions)
#' have been sucessfully called.
#' Accepted value can be the object return by any step function or be feed by
#' \code{\%>\%} from last step function.
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
#' when "pwmfile" is set for motifRc, use this argument to provide
#'  PWM file directory.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome, such as "hg19", "mm10", etc.
#' Default: NULL (e.g. after \code{library (enrichTF)},
#' you can call function \code{setGenome("hg19")})
#' @param  threads \code{Integer} scalar.
#' The maximum threads that will be used in this step.
#' Default: getThreads()
#' @param ... Additional arguments, currently unused.
#' @details
#' Scan for motif occurrences using the prepared PWMs and
#' obtain the promising candidate motifs in these regions.
#' @return An invisible \code{\link{EnrichStep-class}} object
#' (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{genBackground}}
#' \code{\link{findMotifsInRegions}}
#' \code{\link{tfsEnrichInRegions}}
#' @examples
#' setGenome("testgenome") #Use "hg19","hg38",etc. for your application
#' foregroundBedPath <- system.file(package = "enrichTF",
#'     "extdata","testregion.bed")
#' gen <- genBackground(inputForegroundBed = foregroundBedPath)
#' findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate")





setGeneric("enrichFindMotifsInRegions",
           function(prevStep,
                    inputRegionBed = NULL,
                    outputRegionMotifBed = NULL,
                    motifRc = c("integrate","jaspar","pwmfile"),
                    inputPwmFile = getRefFiles("motifpwm"),
                    genome = getGenome(),
                    threads = getThreads(),
                    ...) standardGeneric("enrichFindMotifsInRegions"))



#' @rdname MotifsInRegions
#' @aliases enrichFindMotifsInRegions
#' @export
setMethod(
    f = "enrichFindMotifsInRegions",
    signature = "Step",
    definition = function(prevStep,
                          inputRegionBed = NULL,
                          outputRegionMotifBed = NULL,
                          motifRc = c("integrate","jaspar","pwmfile"),
                          inputPwmFile = getRefFiles("motifpwm"),
                          genome = getGenome(),
                          threads = getThreads(),
                          ...){
        allpara <- c(list(Class = "FindMotifsInRegions",
                          prevSteps = list(prevStep)),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname MotifsInRegions
#' @aliases motifsInRegions
#' @export
findMotifsInRegions <- function(inputRegionBed,
                                outputRegionMotifBed = NULL,
                                motifRc = c("integrate","jaspar","pwmfile"),
                                inputPwmFile = getRefFiles("motifpwm"),
                                genome = getGenome(),
                                threads = getThreads(),
                                ...){
    allpara <- c(list(Class = "FindMotifsInRegions",
                      prevSteps = list()),
                 as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
