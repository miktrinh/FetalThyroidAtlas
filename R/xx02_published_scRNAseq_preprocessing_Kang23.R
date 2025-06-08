## Process published scRNA-seq data of normal adult thyroid tissues

##--- Preprocess the Pu etal 2021 single-cell adult thyroid cancer dataset ---##

outDir = "~/lustre_mt22/Thyroid/Results_v2/xx02_published_scRNAseq_preprocessing"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")

##----------------------------##
##   Set Global parameters  ####
##----------------------------##
maxMT = 30
minGenes = 300
minUMIs = 500  
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = T
scrubScoreMax = 0.5
scrubPath='../cleanCounts/scrubletScores.tsv'
scPath="../cleanCounts/strainedCounts"
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T
rho_max_limit = NULL


##----------------------------------##
##   Preprocessing snRNAseq data  ####
##----------------------------------##
### 1. Import cellranger output data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

outDir_sub = file.path(outDir,'apr25') 

if(is.null(rho_max_limit)){
  plotDir = paste0(file.path(outDir_sub,'Kang23_rhoLimNone_'))
  outPath = paste0(file.path(outDir_sub,'Kang23_rhoLimNone_'))
}

cleanSrat_fp = ifelse(keepMTCells,paste0(outPath,'_clean_withMTCells.RDS'),paste0(outPath,'_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)  
}else{
  if(!dir.exists(outDir_sub)){
    message(sprintf('Creating output directory'))
    dir.create(outDir_sub,recursive = T)
  }
  
  setwd(outDir_sub)
  
  
  
  #### =================================================== #
  srat = NormalizeData(srat,verbose=FALSE)
  # Get cell cycle scores
  data('cc.genes.updated.2019',package='Seurat')
  sPhaseGenes = cc.genes.updated.2019$s.genes
  g2mPhaseGenes = cc.genes.updated.2019$g2m.genes
  #The seed parameter magic is needed because the Seurat authors are dicks
  srat = CellCycleScoring(srat, s.features=sPhaseGenes,g2m.features=g2mPhaseGenes,seed=sample(1e9,1))
  
  # Get MT content
  if(!is.null(mtPattern)){
    srat[['percent.mt']] = PercentageFeatureSet(srat, pattern = mtPattern)
  }
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,rho_max_limit=rho_max_limit,
                      outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)
  
  cleanSrat = QC.output[[1]]
  df.out = QC.output[[2]]
  
  write.csv(df.out,paste0(outPath,'qc_summary.csv'))
}



