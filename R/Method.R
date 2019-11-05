#' @name Pipelines
#' @rdname Pipelines
#' @title ready-to-use pipelines
#' @description There are two ready-to-use pipelines in this package.
#' One is the pipeline mainly for TF enrichment with PECA model.
#' The ohter is the pipeline is for gene regulation analysis including the first pipeline
#' and other related analysis like openness and gene ontology analysis.
#' @param inputForegroundBed \code{Character} scalar.
#' Foreground BED file directory.
#' @param genome \code{Character} scalar.
#' Bioconductor supported genome like "hg19", "mm10", etc.
#' @param threads \code{Numeric} scalar.
#' The max number of threads that can be used by each step of the pipeline
#' @param pipeName \code{Character} scalar or vector.
#' Pipeline name.
#' @param ... Additional arguments to set arguments for each Steps.
#' See below for details.
#' @details This is a function for the pipeline.
#' There are four steps in this pipeline:
#' GenBackground, RegionConnectTarget,
#' FindMotifsInRegions and TFsEnrichInRegions.
#' Parameter setting is available for all these functions. For example,
#' if you want to change
#' the number of background regions (\code{sampleNumb}) into 1000,
#' you can add the argument \code{GenBackground.sampleNumb = 1000}
#' into the function like this:
#' \code{PECA_TF_enrich(inputForegroundBed = "your_file.bed",
#' genome="hg19",GenBackground.sampleNumb = 1000)}.
#' The number of arguments is not limited so you can add other
#' arguments with the format (\code{StepName.argumentName}) in the same way.
#' @author Zheng Wei
#' @return An invisible \code{list} scalar. A list containing all objects that belongs to the pipeline.
#' @references Zhana Duren, et al., Modeling gene regulation from paired
#' expression and chromatin accessibility data.
#' Proc Natl Acad Sci U S A. 2017 1;111(44):15675-80
#' @examples
#'
#' foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
#' # This is the whole pipeline example.
#' PECA_TF_enrich(inputForegroundBed = foregroundBedPath, genome = "testgenome")
#'
#' @export PECA_TF_enrich
PECA_TF_enrich <- function(inputForegroundBed, genome, threads = 2,pipeName = "pipe", ...){
    setGenome(genome)
    setThreads(threads)
    uzp <- unzipAndMergeBed(bedInput = inputForegroundBed, ...)
    gen <- enrichGenBackground(uzp,...)
    conTG <- enrichRegionConnectTargetGene(gen,...)
    findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate",...)
    result <- enrichTFsEnrichInRegions(findMotif,...)
    invisible(getObjsInPipe(pipeName))
}



#' @importFrom magrittr '%>%'
#' @importFrom pipeFrame getObjsInPipe
#' @importFrom pipeFrame setPipeName



#' @aliases Pipelines
#' @rdname Pipelines
#' @export GeneReguPipe
GeneReguPipe <- function(inputForegroundBed, genome, threads = 2, pipeName = "pipe", ...){
    setGenome(genome)
    setThreads(threads)
    setPipeName(pipeName = pipeName)
    obj <- unzipAndMergeBed(bedInput = inputForegroundBed, ...) %>%
        enrichTissueOpennessConserve(...)%>%
        enrichTissueOpennessSpecificity(...)%>%
        enrichGenBackground(...) %>%
        enrichRegionConnectTargetGene(...) %>%
        enrichFindMotifsInRegions(motifRc="integrate",...) %>%
        enrichTFsEnrichInRegions(...) %>%
        enrichGeneOntology(...)%>%
        enrichSingleSampleReport(...)
    invisible(getObjsInPipe(pipeName))
}



