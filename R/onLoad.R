
#' @importFrom pipeFrame Step
#' @importFrom pipeFrame addEdges
#' @importFrom pipeFrame checkAndInstallBSgenome
#' @importFrom pipeFrame checkAndInstallGenomeFa
#' @importFrom pipeFrame checkFileCreatable
#' @importFrom pipeFrame checkFileExist
#' @importFrom pipeFrame checkPathExist
#' @importFrom pipeFrame getBasenamePrefix
#' @importFrom pipeFrame getGenome
#' @importFrom pipeFrame getJobDir
#' @importFrom pipeFrame getJobName
#' @importFrom pipeFrame getNextSteps
#' @importFrom pipeFrame getPathPrefix
#' @importFrom pipeFrame getPrevSteps
#' @importFrom pipeFrame getRef
#' @importFrom pipeFrame getRefDir
#' @importFrom pipeFrame getRefFiles
#' @importFrom pipeFrame getRefRc
#' @importFrom pipeFrame getThreads
#' @importFrom pipeFrame getTmpDir
#' @importFrom pipeFrame getValidGenome
#' @importFrom pipeFrame initPipeFrame
#' @importFrom pipeFrame processing
#' @importFrom pipeFrame runWithFinishCheck
#' @importFrom pipeFrame setGenome
#' @importFrom pipeFrame setJobName
#' @importFrom pipeFrame setRefDir
#' @importFrom pipeFrame setThreads
#' @importFrom pipeFrame setTmpDir
#' @importMethodsFrom pipeFrame checkAllPath
#' @importMethodsFrom pipeFrame checkRequireParam
#' @importMethodsFrom pipeFrame clearStepCache
#' @importMethodsFrom pipeFrame getAutoPath
#' @importMethodsFrom pipeFrame stepName
#' @importMethodsFrom pipeFrame getParam
#' @importMethodsFrom pipeFrame getParamItems
#' @importMethodsFrom pipeFrame getParamMD5Path
#' @importMethodsFrom pipeFrame stepID
#' @importMethodsFrom pipeFrame getStepWorkDir
#' @importMethodsFrom pipeFrame init
#' @importMethodsFrom pipeFrame isReady
#' @importMethodsFrom pipeFrame writeLog
.onLoad <- function(libname, pkgname) {
    initPipeFrame(availableGenome = c("hg19", "hg38", "mm9",
                                      "mm10","testgenome"),
                  defaultJobName = paste0(pkgname,"-pipeline"),
                  defaultCheckAndInstallFunc = checkAndInstall
    )

    addEdges(edges = c( "UnzipAndMergeBed","GenBackground",
                        "UnzipAndMergeBed","TissueOpennessConserve",
                        "UnzipAndMergeBed","TissueOpennessSpecificity"),
             argOrder = 1)
    addEdges(edges = c( "FindMotifsInRegions","TFsEnrichInRegions"),
             argOrder = 2)
    addEdges(edges = c("RegionConnectTargetGene","TFsEnrichInRegions"),
             argOrder = 3)
    addEdges(edges = c("GenBackground","FindMotifsInRegions",
                       "GenBackground","RegionConnectTargetGene",
                       "GenBackground","TFsEnrichInRegions"),
             argOrder = 1)



}
