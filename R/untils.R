
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

dowloadChromSize <- function(resultDirPaths){
    download.file(url = sprintf("https://wzthu.github.io/esATAC/refdata/%s.DHS.bed.gz",getGenome()),
                  destfile = paste0(resultDirPaths,".gz"),method = getOption("download.file.method"))
    gunzip(paste0(resultDirPaths,".gz"),destname=resultDirPaths,overwrite=TRUE)

}
checkAndInstall <- function(check = TRUE, ...){
    runWithFinishCheck(func = checkAndInstallBSgenome,refName = "bsgenome", resultVal = getBSgenome(getGenome()), execForNonRsFile = check)
    runWithFinishCheck(func = checkAndInstallGenomeFa,refName = "fasta", resultDirPaths = paste0(getGenome(),".fa"))
    runWithFinishCheck(func = dowloadChromSize,refName = "DHS", resultDirPaths = "DHS.bed")
}


