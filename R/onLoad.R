.onLoad <- function(libname, pkgname) {
    initPipeFrame(availableGenome = c("hg19", "hg38", "mm9", "mm10"),
                  defaultJobName = paste0(pkgname,"-pipeline"),
                  defaultCheckAndInstallFunc = checkAndInstall
    )


    addEdges(edges = c( "FindMotifsInRegions","TFsEnrichInRegions"),
             argOrder = 2)
    addEdges(edges = c("RegionConnectTargetGene","TFsEnrichInRegions"),
             argOrder = 3)
    addEdges(edges = c("GenBackground","FindMotifsInRegions",
                       "GenBackground","RegionConnectTargetGene",
                       "GenBackground","TFsEnrichInRegions"),
             argOrder = 1)



}
