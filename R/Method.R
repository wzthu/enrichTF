PECA_TF_enrich <- function(bedInput,genome,...){
    setGenome(genome)
    gen <- genBackground(inputForegroundBed = bedInput,...)
    conTG <- enrichRegionConnectTargetGene(gen,...)

    findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate",...)

    result <- enrichTFsEnrichInRegions(gen,findMotif,conTG,...)
    return(list(gen,conTG,findMotif,result))
}
