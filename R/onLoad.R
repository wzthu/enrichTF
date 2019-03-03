.onLoad <- function(libname, pkgname) {
    initPipeFrame(availableGenome = c("hg19", "hg38", "mm9", "mm10"),
                  defaultJobName = paste0(pkgname,"-pipeline"),
                  defaultCheckAndInstallFunc = checkAndInstall
    )
}
