
checkAndInstallBSgenome <- function(resultDirPaths){
    genome <- getGenome()
    bsgenomename<- BSgenome::available.genomes()[grepl(paste0(genome,"$"),BSgenome::available.genomes())]
    if(length(bsgenomename)==0){
        message()
        stop("There is no BSgenome support for this genome")
    }
    bsgenomeinstall <- BSgenome::installed.genomes()[grepl(paste0(genome,"$"),BSgenome::installed.genomes())]
    if(length(bsgenomeinstall)==0){
        message(paste("BSgenome for ",genome,"has not been installed,"))
        message("begin to install ...")
        BiocManager::install(bsgenomename)
    }
}

checkAndInstallGenomeFa <- function(resultDirPaths){
    outFile <- resultDirPaths
    bsgenome<-getRefRs("bsgenome")
    if(!is(bsgenome, "BSgenome")){stop("The variable 'bsgenome' is not a BSgenome")}
    append <- FALSE
    for(chrT in seqnames(bsgenome)){
        if(is.null(masks(bsgenome[[chrT]])))
            chrSeq <- DNAStringSet(bsgenome[[chrT]])
        else
            chrSeq <- DNAStringSet(injectHardMask(bsgenome[[chrT]], letter="N"))
        names(chrSeq) <- chrT
        writeXStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
        append <- TRUE
    }
    return(outFile)
}

dowloadMotifFile <- function(resultDirPaths){
    download.file(url = "https://wzthu.github.io/enrich/refdata/all_motif_rmdup",
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}

dowloadREgeneFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = sprintf("https://wzthu.github.io/enrich/refdata/%s/RE_gene_corr_hg19.bed",genome),
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}


dowloadEnhancerREgeneFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = sprintf("https://wzthu.github.io/enrich/refdata/%s/Enhancer_RE_gene_corr_hg19.bed",genome),
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}

convertPWMFileToPWMobj <- function(resultDirPaths){
    motiffile <-getRefFiles("motifpwm")
    pwm <- getMotifInfo1(motiffile)
    save(pwm,file = resultDirPaths)
}


dowloadMotifTFTableFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = "https://wzthu.github.io/enrich/refdata/MotifTFTable.RData",
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}
dowloadMotifWeightsFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = "https://wzthu.github.io/enrich/refdata/MotifWeights.RData",
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}
dowloadTFgeneRelMtxFile <- function(resultDirPaths){
    genome <-getGenome()
    download.file(url = "https://wzthu.github.io/enrich/refdata/TFgeneRelMtx.RData",
                  destfile = resultDirPaths,method = getOption("download.file.method"))

}

checkAndInstall <- function(check = TRUE, ...){
    runWithFinishCheck(func = checkAndInstallBSgenome,refName = "bsgenome", resultVal = getBSgenome(getGenome()), execForNonRsFile = check)
    runWithFinishCheck(func = checkAndInstallGenomeFa,refName = "fasta", resultDirPaths = paste0(getGenome(),".fa"))
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
                                                ID = motifSrc,
                                                tags = list(seq=exseq,
                                                            motifName = motifName,
                                                            motifSrc=motifSrc,
                                                            motifPlf = motifPlf,
                                                            threshold = as.numeric(threshold)))
                PWMList[[motifName]] <- p_matrix

            }
            strlist <- unlist(unlist(strsplit(substring(line,2),"\t|/")))
            exseq <- strlist[1]
            motifName <- strlist[2]
            motifSrc <- strlist[3]
            motifPlf <- strlist[4]
            threshold <- strlist[5]
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
                                    ID = motifSrc,
                                    tags = list(seq=exseq,
                                                motifName = motifName,
                                                motifSrc=motifSrc,
                                                motifPlf = motifPlf,
                                                threshold = as.numeric(threshold)))
    PWMList[[motifName]] <- p_matrix

    PWMList <- do.call(TFBSTools::PWMatrixList, PWMList)
    print(names(PWMList))
    return(PWMList)
}
