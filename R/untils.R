#' @importFrom R.utils gunzip
#' @importFrom pipeFrame getRefRc
#' @importFrom pipeFrame checkAndInstallBSgenome
#' @importFrom pipeFrame checkAndInstallGenomeFa
#' @importFrom pipeFrame getRefFiles
#' @importFrom BSgenome getBSgenome


dowloadMotifFile <- function(resultDirPaths){
    download.file(url = "https://wzthu.github.io/enrich/refdata/all_motif_rmdup",
                  destfile = resultDirPaths)

}

dowloadREgeneFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = sprintf("https://wzthu.github.io/enrich/refdata/%s/RE_gene_corr_hg19.bed",genome),
                  destfile = resultDirPaths)

}


dowloadEnhancerREgeneFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = sprintf("https://wzthu.github.io/enrich/refdata/%s/Enhancer_RE_gene_corr_hg19.bed",genome),
                  destfile = resultDirPaths)

}

convertPWMFileToPWMobj <- function(resultDirPaths){
    motiffile <-getRefFiles("motifpwm")
    pwm <- getMotifInfo1(motiffile)
    save(pwm,file = resultDirPaths)
}


dowloadMotifTFTableFile <- function(resultDirPaths){
    download.file(url = "https://wzthu.github.io/enrich/refdata/MotifTFTable.RData",
                  destfile = resultDirPaths)

}
dowloadMotifWeightsFile <- function(resultDirPaths){
    download.file(url = "https://wzthu.github.io/enrich/refdata/MotifWeights.RData",
                  destfile = resultDirPaths)

}
dowloadTFgeneRelMtxFile <- function(resultDirPaths){
    download.file(url = "https://wzthu.github.io/enrich/refdata/TFgeneRelMtx.RData",
                  destfile = resultDirPaths)

}

checkAndInstall <- function(check = TRUE, ...){
    runWithFinishCheck(func = checkAndInstallBSgenome,refName = "bsgenome", resultVal = getBSgenome(getGenome()), execForNonRsFile = check)
#    runWithFinishCheck(func = checkAndInstallGenomeFa,refName = "fasta", resultDirPaths = paste0(getGenome(),".fa"))
    runWithFinishCheck(func = dowloadMotifFile,refName = "motifpwm", resultDirPaths = "motifpwm")
    runWithFinishCheck(func = convertPWMFileToPWMobj, "motifPWMOBJ", resultDirPaths = "motifPWMOBJ.RData")
    runWithFinishCheck(func = dowloadREgeneFile, "RE_gene_corr", resultDirPaths = "RE_gene_corr.bed")
    runWithFinishCheck(func = dowloadEnhancerREgeneFile, "Enhancer_RE_gene_corr", resultDirPaths = "Enhancer_RE_gene_corr.bed")
    runWithFinishCheck(func = dowloadMotifTFTableFile, "MotifTFTable", resultDirPaths = "MotifTFTable.RData")
    runWithFinishCheck(func = dowloadMotifWeightsFile, "MotifWeights", resultDirPaths = "MotifWeights.RData")
    runWithFinishCheck(func = dowloadTFgeneRelMtxFile, "TFgeneRelMtx", resultDirPaths = "TFgeneRelMtx.RData")
}


getMotifInfo1 <- function(motiffile = NULL){
    PWMList <- list()


    con <- file(motiffile, "r")
    lines <- readLines(con)
    close(con)
    exseq <- NULL
    motifName <- NULL
    motifSrc <- NULL
    motifPlf <- NULL
    threshold <- NULL
    p <- c()
    for(line in lines){
        if(substring(line, 1, 1) == ">"){
            if(!is.null(motifName)){
                pwm <- matrix(data = p, nrow = 4, dimnames = list(c("A","C","G","T")))
                p_matrix <- TFBSTools::PWMatrix(profileMatrix = pwm,
                                                name = motifName,
                                                tags = list(seq=exseq,
                                                            motifName = motifName,
                                                            threshold = as.numeric(threshold)))
                PWMList[[motifName]] <- p_matrix

            }
            strlist <- unlist(unlist(strsplit(substring(line,2),"\t")))
            exseq <- strlist[1]
            motifName <- strlist[2]
            threshold <- strlist[3]
            p <- c()
        }else{
            val <- as.numeric(unlist(strsplit(line,"\t")))
            val <- val/sum(val)
            p <- c(p,val)
        }
    }
    pwm <- matrix(data = as.numeric(p), nrow = 4, dimnames = list(c("A","C","G","T")))
    p_matrix <- TFBSTools::PWMatrix(profileMatrix = pwm,
                                    name = motifName,
                                    tags = list(seq=exseq,
                                                motifName = motifName,
                                                threshold = as.numeric(threshold)))
    PWMList[[motifName]] <- p_matrix

    PWMList <- do.call(TFBSTools::PWMatrixList, PWMList)
    print(names(PWMList))
    return(PWMList)
}
