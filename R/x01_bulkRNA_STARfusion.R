## Process STAR-fusion output to check for RET fusion in each bulk sample

outDir = '~/lustre_mt22/Thyroid/Results_v2/xx04_bulkRNA_STARfusion'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)




##-------------##
##  Libraries  ##
##-------------##
library(tidyverse)


##---- paediatric bulk RNAseq
sce_path = '~/lustre_mt22/Thyroid/Data/inhouse_bulkRNA_thyroid/inhouse_bulkRNA_thyroid_2410_sce.RDS'
inhouse_sce = readRDS(sce_path)
inhouse_mdat = as.data.frame(colData(inhouse_sce))
inhouse_mdat = inhouse_mdat[match(colnames(inhouse_rawCnt),inhouse_mdat$sampleID),]



STAR_fusion_outDir = '/lustre/scratch127/cellgen/cellgeni/tickets/tic-3426/results/STAR-Fusion/'
STAR_fusion_outFiles = list.files(STAR_fusion_outDir,recursive = T,pattern = 'star-fusion.fusion_predictions.tsv',full.names = T)
names(STAR_fusion_outFiles) = gsub('_fusion','',basename(dirname(dirname(STAR_fusion_outFiles))))

allFusion = data.frame()
for(i in 1:length(STAR_fusion_outFiles)){
  sample = names(STAR_fusion_outFiles)[i]
  fusion_result = read.delim(STAR_fusion_outFiles[i],sep = '\t')
  if(nrow(fusion_result) == 0){
    print(inhouse_mdat$age[inhouse_mdat$sangerSampleID == sample])
    next
  }
  fusion_result$sampleID = sample
  if(nrow(allFusion) == 0 ){
    allFusion = fusion_result
  }else{
    allFusion = rbind(allFusion,fusion_result)  
  }
}

## Add sample level metadata
allFusion = cbind(allFusion,inhouse_mdat[match(allFusion$sampleID,inhouse_mdat$sangerSampleID),!colnames(inhouse_mdat) %in% colnames(allFusion)])
write.csv(allFusion,'bulk_paed_foetal_thyroid_STAR_fusion_output_allSamples_combined.csv')
allFusion = read.csv('bulk_paed_foetal_thyroid_STAR_fusion_output_allSamples_combined.csv')
View(inhouse_mdat[inhouse_mdat$sangerSampleID %in% allFusion$sampleID[grepl('RET',allFusion$X.FusionName)],])

