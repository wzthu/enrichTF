#' @name PECA_TF_enrich
#' @title TF enrichment with PECA model
#' @description This is a pipeline for TF enrichment with PECA model.
#' @param inputForegroundBed \code{Character} scalar.
#' Foreground BED file directory.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome like "hg19", "mm10", etc.
#' @param threads \code{Numeric} scalar.
#' The max number of threads that can be used by each step of the pipeline
#' @param ... Additional arguments to set arguments for each Steps.
#' See below for details.
#' @details This is a function for the pipeline. There are four steps in this pipeline: GenBackground, RegionConnectTarget,
#' FindMotifsInRegions and TFsEnrichInRegions. Parameter setting is available for all these functions. For example, if you want to change
#' the number of background regions (\code{sampleNumb}) into 1000,
#' you can add the argument \code{GenBackground.sampleNumb = 1000} into the function like this:
#' \code{PECA_TF_enrich(inputForegroundBed = "your_file.bed",genome="hg19",GenBackground.sampleNumb = 1000)}.
#' The number of arguments is not limited so you can add other arguments with the format (\code{StepName.argumentName}) in the same way.
#' @author Zheng Wei
#' @return An invisible \code{list} contains all four steps \code{EnrichTF} objects
#' @references Zhana Duren, et al., Modeling gene regulation from paired expression and chromatin accessibility data.
#' Proc Natl Acad Sci U S A. 2017 1;111(44):15675-80
#' @export PECA_TF_enrich
#' @examples
#'
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' # This is the whole pipeline example.
#' # PECA_TF_enrich(inputForegroundBed = foregroundBedPath, genome = "testgenome")
#'

PECA_TF_enrich <- function(inputForegroundBed, genome, threads = 2, ...){
    setGenome(genome)
    setThreads(threads)
    gen <- genBackground(inputForegroundBed = inputForegroundBed,...)
    conTG <- enrichRegionConnectTargetGene(gen,...)
    findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate",...)
    result <- enrichTFsEnrichInRegions(gen,findMotif,conTG,...)
    invisible(list(GenBackground = gen,
                   RegionConnectTargetGene = conTG,
                   FindMotifsInRegions = findMotif,
                   TFsEnrichInRegions = result))
}
