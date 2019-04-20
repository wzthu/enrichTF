#' @importFrom R.utils gunzip
#' @importFrom pipeFrame getRefRc
#' @importFrom pipeFrame checkAndInstallBSgenome
#' @importFrom pipeFrame checkAndInstallGenomeFa
#' @importFrom pipeFrame getRefFiles
#' @importFrom BSgenome getBSgenome
#' @importFrom utils download.file
#' @importFrom R.utils gunzip
#' @importFrom pipeFrame input "input<-" output "output<-" "param" "param<-"


downloadAndGunzip <- function(urlplaceholder,refFilePath){
    genome <- getGenome()
    download.file(url = sprintf(urlplaceholder,genome),
                  destfile = paste0(refFilePath,'.gz'))
    gunzip(paste0(refFilePath,'.gz'),remove = TRUE)
}

copyAndGunzip <- function(fileName,refFilePath){
    file.copy(from = system.file(package = "enrichTF", "extdata",fileName),
                    to = paste0(refFilePath,'.gz'),overwrite = TRUE)
    gunzip(paste0(refFilePath,'.gz'),remove = TRUE)
}

dowloadMotifFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("all_motif_rmdup.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/all_motif_rmdup.gz",
            refFilePath)
    }
}

dowloadREgeneFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("RE_gene_corr.bed.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/RE_gene_corr.bed.gz",
            refFilePath)
    }
}


dowloadEnhancerREgeneFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("Enhancer_RE_gene_corr.bed.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/Enhancer_RE_gene_corr.bed.gz",
            refFilePath)
    }

}

convertPWMFileToPWMobj <- function(refFilePath){
    motiffile <-getRefFiles("motifpwm")
    pwm <- getMotifInfo1(motiffile)
    save(pwm,file = refFilePath)
}


dowloadMotifTFTableFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("MotifTFTable.RData.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/MotifTFTable.RData.gz",
            refFilePath)
    }
}
dowloadMotifWeightsFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("MotifWeights.RData.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/MotifWeights.RData.gz",
                          refFilePath)
    }
}
dowloadTFgeneRelMtxFile <- function(refFilePath){
    if(getGenome() == "testgenome"){
        copyAndGunzip("TFgeneRelMtx.RData.gz",refFilePath)
    }else{
        downloadAndGunzip(
            "https://wzthu.github.io/enrich/refdata/%s/TFgeneRelMtx.RData.gz",
            refFilePath)
    }
}

checkAndInstallBSgenomeTestgenome <- function(refFilePath){
    genome <- getGenome()
    if(genome == "testgenome"){
        genome <- "hg19"
    }
    checkAndInstallBSgenome(refFilePath, genome)
}

checkAndInstall <- function(check = TRUE, ...){
    runWithFinishCheck(func = checkAndInstallBSgenomeTestgenome,
                       refName = "bsgenome")
#    runWithFinishCheck(func = checkAndInstallGenomeFa,refName = "fasta", refFilePath = paste0(getGenome(),".fa"))
    runWithFinishCheck(func = dowloadMotifFile,
                       refName = "motifpwm",
                       refFilePath = "motifpwm")
    runWithFinishCheck(func = convertPWMFileToPWMobj,
                       "motifPWMOBJ",
                       refFilePath = "motifPWMOBJ.RData")
    runWithFinishCheck(func = dowloadREgeneFile,
                       "RE_gene_corr",
                       refFilePath = "RE_gene_corr.bed")
    runWithFinishCheck(func = dowloadEnhancerREgeneFile,
                       "Enhancer_RE_gene_corr",
                       refFilePath = "Enhancer_RE_gene_corr.bed")
    runWithFinishCheck(func = dowloadMotifTFTableFile,
                       "MotifTFTable",
                       refFilePath = "MotifTFTable.RData")
    runWithFinishCheck(func = dowloadMotifWeightsFile,
                       "MotifWeights",
                       refFilePath = "MotifWeights.RData")
    runWithFinishCheck(func = dowloadTFgeneRelMtxFile,
                       "TFgeneRelMtx",
                       refFilePath = "TFgeneRelMtx.RData")
}


getMotifInfo1 <- function(motiffile = NULL){
    #PWMList <- list()


    con <- file(motiffile, "r")
    lines <- readLines(con)
    close(con)
    exseq <- NULL
    motifName <- NULL
    motifPlf <- NULL
    threshold <- NULL
    p <- c()
    # for(line in lines){
    #     if(substring(line, 1, 1) == ">"){
    #         if(!is.null(motifName)){
    #             pwm <- log2(matrix(data = p, nrow = 4, dimnames = list(c("A","C","G","T"))) * 4)
    #             p_matrix <- TFBSTools::PWMatrix(profileMatrix = pwm,
    #                                             name = motifName,
    #                                             tags = list(seq=exseq,
    #                                                         motifName = motifName,
    #                                                         threshold = as.numeric(threshold)))
    #             PWMList[[motifName]] <- p_matrix
    #
    #         }
    #         strlist <- unlist(unlist(strsplit(substring(line,2),"\t")))
    #         exseq <- strlist[1]
    #         motifName <- strlist[2]
    #         threshold <- strlist[3]
    #         p <- c()
    #     }else{
    #         val <- as.numeric(unlist(strsplit(line,"\t")))
    #         val <- val/sum(val)
    #         p <- c(p,val)
    #     }
    # }
    rs <- lapply(lines, function(line){
        if(substring(line, 1, 1) == ">"){
            strlist <- unlist(unlist(strsplit(substring(line,2),"\t")))
            exseq <- strlist[1]
            motifName <- strlist[2]
            threshold <- strlist[3]
            return(list(p=NULL, exseq = exseq,
                        motifName = motifName, threshold = threshold))
        }else{
            val <- as.numeric(unlist(strsplit(line,"\t")))
            val <- val/sum(val)
            return(list(p=val, exseq = NULL,
                        motifName = NULL, threshold = NULL))
        }
    })
    cutpoint <- lapply(seq_len(length(rs)), function(s){
        if(is.null(rs[[s]][["motifName"]])){
           return(NULL)
        }else{return(s)}
    })
    cutpoint <- unlist(cutpoint)

    cutpoint <- c(cutpoint,length(lines) + 1 )
    PWMList <- lapply(seq_len(length(cutpoint)-1), function(i){
        begin <- cutpoint[i]
        exseq <- rs[[begin]][["exseq"]]
        motifName <- rs[[begin]][["motifName"]]
        threshold <- rs[[begin]][["threshold"]]
        begin <- begin + 1
        end <- cutpoint[i + 1] - 1
        p <- lapply(rs[begin:end], function(t){
            return(t$p)
        })
        p <- do.call("c",p)
        pwm <- log2(matrix(data = p, nrow = 4, dimnames = list(c("A","C","G","T"))) * 4)
        p_matrix <- TFBSTools::PWMatrix(profileMatrix = pwm,
                                            name = motifName,
                                            tags = list(seq=exseq,
                                                        motifName = motifName,
                                                        threshold = as.numeric(threshold)))
        pwmrs <- list()
        pwmrs[[motifName]] <- p_matrix
        return(p_matrix)
    })
    PWMList <- do.call("c",PWMList)
    PWMList <- do.call(TFBSTools::PWMatrixList, PWMList)
    return(PWMList)
}
