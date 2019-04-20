context("test whole pipeline")


test_that("test whole pipeline",{


    library(magrittr)
    setGenome("testgenome") #Use "hg19","hg38",etc. for your application
    foregroundBedPath <- system.file(package = "enrichTF", "extdata","testregion.bed")
    gen <- genBackground(inputForegroundBed = foregroundBedPath)
    conTG <- enrichRegionConnectTargetGene(gen)
    findMotif <- enrichFindMotifsInRegions(gen,motifRc="integrate")
    result <- enrichTFsEnrichInRegions(gen)

    expect_true(file.exists(gen$outputForegroundBed))
    expect_true(file.exists(gen$outputBackgroundBed))
    expect_true(file.exists(gen$outputRegionBed))
    expect_true(file.exists(conTG$outputForegroundBed))
    expect_true(file.exists(conTG$outputBackgroundBed))
    expect_true(file.exists(conTG$regularGeneCorrBed))
    expect_true(file.exists(conTG$enhancerRegularGeneCorrBed))
    expect_true(file.exists(findMotif$inputRegionBed))
    expect_true(file.exists(findMotif$outputRegionMotifBed))
    expect_true(file.exists(result$outputTFsEnrichTxt))



})
