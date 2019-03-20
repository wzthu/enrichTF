#' @name PECA_TF_enrich
#' @title TF enrich with PECA2 mothod
#' @description TF enrich with PECA2 mothod
#' @param inputForegroundBed \code{Character} scalar.
#' Foreground BED file directory.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome like "hg19", "mm10", etc.
#' @param threads \code{Numeric} scalar.
#' The max number of threads that can be used by each step of the pipeline
#' @param ... Additional arguments to set arguments for each Steps.
#' See below for details.
#' @details There four steps in this pipeline: GenBackground, RegionConnectTarget,
#' FindMotifsInRegions and TFsEnrichInRegions. For instance, if you want to change
#' the number of background regions (\code{sampleNumb}) into 1000,
#' you can add argument \code{GenBackground.sampleNumb = 1000} to the function in this way:
#' \code{PECA_TF_enrich(inputForegroundBed = "your/file.bed",genome="hg19",GenBackground.sampleNumb = 1000)}.
#' The number of arguments is not limited so you can add other arguments with \code{StepName.argumentName} in the same way.
#' @author Zheng Wei
#' @return An invisible \code{list} contain all four steps \code{EnrichTF} objects
#' @references Duren Z.,.....PNAS
#' @export PECA_TF_enrich
#' @examples
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' PECA_TF_enrich(inputForegroundBed = foregroundBedPath, genome = "testgenome")

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
