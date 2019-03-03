#
# randomSampleOnGenome(rangeLen, sampleNumber){
#     bsgenome <- getRefRs("bsgenome")
#     chrlens <-seqlengths(bsgenome)
#     selchr <- grep("_|M",names(chrlens),invert=TRUE)
#     chrlens <- chrlens[selchr]
#     startchrlens <- chrlens - rangeLen
#     totallen <- sum(startchrlens)
#     spnb <- floor(runif(sampleNumber) * totallen) + 1
#     acclen <- startchrlens
#     for(i in 2:length(acclen)){
#         acclen[i] <- acclen[i-1] + acclen[i]
#     }
#     acclen <- c(0,acclen)
#     gr <- GRanges()
#     for(i in 1:(length(acclen)-1)){
#         sel <- spnb[(spnb>acclen[i]) & (spnb<=acclen[i+1])]
#         sel <- sel - acclen[i]
#         gr <- c(gr,GRanges(seqnames = names(acclen)[i+1], ranges = IRanges(start = sel, width = 1000)))
#     }
#     return(sort(gr,ignore.strand=TRUE))
# }
#
#
# genBackground <- function(inputRegion, outputForeground = NULL, outputBackgroudPath = NULL, outputRegionBedPath = NULL, rangeLen =1000, backgroundnb = 10000){
#     foregroundgr <- import(con = inputRegionBedPath, genome = getGenome(),format = "bed")
#     midpoint <- (start(inputBed) + end(inputBed))/2
#     start(foregroundgr) <- floor(midpoint - rangeLen/2)
#     end(foregroundgr) <- floor(midpoint + rangeLen/2)
#     foregroundgr <- sort(foregroundgr,ignore.strand=TRUE)
#     backgroundgr <- randomSampleOnGenome(rangeLen, backgroundnb)
#     regiongr <- c(foregroundgr,backgroundgr)
# }
#
#
# motifScan <- function(regions, motifTarget){
#     opts <- list()
#     opts[["tax_group"]] <- "vertebrates"
#     pfm <- getMatrixSet(JASPAR2018::JASPAR2018, opts)
#     names(pfm) <- TFBSTools::name(pfm)
#     #    names(pfm) <- gsub(pattern = "[^a-zA-Z0-9:)(]", replacement = "", x = TFBSTools::name(pfm), perl = TRUE)
#     data(example_motifs, package = "motifmatchr")
#     regions <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#                                       ranges = IRanges::IRanges(start = c(76585873,42772928,
#                                                                           100183786),
#                                                                 width = 500))
#     motif_ix <-matchMotifs(pwms = example_motifs, subject = regions, genome = getGenome(),out="positions")
#
#     gr <- GRanges()
#     for(i in 1:length(motif_ix)){
#         gr<-c(gr,subsetByOverlaps(regions,motif_ix[[i]]))
#     }
#     return(gr)
# }
#
#
# groudGeneOverlap <- function(gr,gen){
#     a<-subsetByOverlaps(gr,gen)
#     mcols(a)<- mcols(a)[1]
#     return(a)
# }
#
# connectTargetGene <- function(forgroundBed, backgroundBed){
#     tgfiles<-getRefFiles("targetGene")
#     regene <- import(con = tgfiles[1], genome = getGenome(),format = "bed" )
#     enhancerregene <- import(con = tgfiles[2], genome = getGenome(),format = "bed" )
#     fr <- groundGeneOverlap(forgroundBed,regene)
#     br <- groundGeneOverlap(backgroundBed,regene)
#     fe <- groundGeneOverlap(forgroundBed,enhancerregene)
#     be <- groundGeneOverlap(backgroundBed,enhancerregene)
#     export.bed(fr)
#     export.bed(br)
#     export.bed(fe)
#     export.bed(be)
# }
#
#
# enrichAnalysis <- function(){
#
# }
#
#
#
#
#
#
# test <- function() {
#     setRefdir()
#     setGenome("hg19")
#     genBackground()
#
# }
