

# Helper functions
import_bulkRNA_thyroid = function(bulk_sources = c('TCGA_Thyroid','inhouse'),inhouse_bulk_fp = '/lustre/scratch126/cellgen/behjati/mt22/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS'){
  library(readxl)
  library(SummarizedExperiment)
  library(edgeR)
  library(singscore)
  library(GenomicFeatures)
  
  
  ##----- 1. Import bulk counts -----##
  bulk_samples = tibble()
  sce_list = list()
  raw_count = data.frame()
  tpmCnt = data.frame()
  
  if('inhouse' %in% bulk_sources){
    se_path = inhouse_bulk_fp
    inhouse_se = readRDS(se_path)
    sce_list[['Sanger']] = inhouse_se
    inhouse_mdat = as.data.frame(colData(inhouse_se))
    
    bulk_samples = rbind(bulk_samples,inhouse_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    inhouse_rawCnt = assays(inhouse_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = inhouse_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(inhouse_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],inhouse_rawCnt[genesToKeep,])
    }
  }
  
  print('Getting TPM counts')
  ## Get TPM count
  tpmCnt = data.frame()
  for(i in 1:length(sce_list)){
    print(i)
    tpm = assays(sce_list[[i]])[['counts_tpm']]
    # if(names(sce_list)[i] == 'GTEx_Thyroid'){
    #   rowDat = rowData(sce_list[[i]])
    #   rownames(tpm) = rowDat$ensID[match(rownames(tpm),rowDat$gene_id)]
    # }
    
    if(nrow(tpmCnt) == 0){
      tpmCnt = tpm
    }else{
      genesToKeep = intersect(rownames(tpmCnt),rownames(tpm))
      print(length(genesToKeep))
      tpmCnt = cbind(tpmCnt[genesToKeep,],tpm[genesToKeep,])
    }
  }
  
  
  
  ##----- 2. Import bulk metadata -----##
  ##---- Fix some spelling / Categorise age ----##
  bulk_samples$ageCat = ifelse(bulk_samples$source == 'fAdrenal_SCPs','0',
                               ifelse(bulk_samples$age %in% c('4.3','4.5','5.1','5.2'),'4-5',
                                      ifelse(bulk_samples$age %in% c('6.4','7.1','7.4','7.6','7.8'),'6-7',
                                             ifelse(bulk_samples$age %in% c('8.9','9','9.6','9.9'),'8-9',
                                                    ifelse(bulk_samples$age %in% c('13.3','13.7','14'),'13-14',
                                                           ifelse(bulk_samples$age %in% as.character(seq(18:30)),'18-30',
                                                                  ifelse(bulk_samples$age %in% as.character(seq(31:60)),'31-60',
                                                                         ifelse(bulk_samples$age %in% as.character(seq(61:81)),'61-81',bulk_samples$age))))))))
  
  
  bulk_samples = bulk_samples[!bulk_samples$cancerType %in% c('PLEUROPULMONARY BLASTOMA'),]
  
  # # Add full path to bulk_samples$sample_file
  # bulk_samples$sample_file_fullpath = paste0('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results_v2/x05_bulkCSA/bulkData/',bulk_samples$sampleID,'_counts.txt')
  
  ## cancerType_group
  bulk_samples$cancerType_details = bulk_samples$cancerType
  bulk_samples$cancerType[grepl('Normal',bulk_samples$cancerType)] = 'Normal'
  bulk_samples$cancerType[grepl('Normal.adj',bulk_samples$cancerType_details)] = 'Normal.adj'
  bulk_samples$cancerType[grepl('ollicular|FOLLICULAR|fvPTC|miFTC',bulk_samples$cancerType)] = 'follicular_PTC'
  bulk_samples$cancerType[grepl('Tall Cell',bulk_samples$cancerType)] = 'tallCell_PTC'
  bulk_samples$cancerType[grepl('other|Other|PLEUROPULMONARY BLASTOMA',bulk_samples$cancerType)] = 'others'
  bulk_samples$cancerType[bulk_samples$cancerType %in% c('CCDC6_RET:Classical/usual')] = 'CCDC6_RET_PTC_TCGA'
  bulk_samples$cancerType[bulk_samples$cancerType %in% c('PTC_BRAF') & bulk_samples$source == 'Lee_2021'] = 'BRAF_PTC_paediatric'
  bulk_samples$cancerType[bulk_samples$cancerType %in% c('PTC_fusion') & bulk_samples$source == 'Lee_2021'] = 'fusion_PTC_paediatric'
  
  bulk_samples$cancerType[grepl('PAPILLARY CARCINOMA|cPTC|Papillary thyroid carcinoma|PTC',bulk_samples$cancerType) & bulk_samples$source %in% c('stJudes_Thyroid','snRNAseq_Y24')] = 'PTC_paediatric'
  bulk_samples$cancerType[grepl('-RET|_RET',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'RETfusion_PTC_paediatric'
  bulk_samples$cancerType[grepl('-BRAF',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'BRAFfusion_PTC_paediatric'
  bulk_samples$cancerType[grepl('-NTRK3',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'NTRK3fusion_PTC_paediatric'
  bulk_samples$cancerType[grepl('-ALK|_ALK',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'ALKfusion_PTC_paediatric'
  bulk_samples$cancerType[grepl(':Not Available',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'PTC_paediatric'
  bulk_samples$cancerType[grepl(':BRAF|:HRAS|DICER1',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'somMut_PTC_paediatric'
  bulk_samples$cancerType[grepl(':None',bulk_samples$cancerType_details) & bulk_samples$cancerType == 'PTC_paediatric' & bulk_samples$source == 'stJudes_Thyroid'] = 'uninformed_PTC_paediatric'
  
  bulk_samples$cancerType[grepl('Primary Tumor|Classical/usual|cPTC',bulk_samples$cancerType)] = 'PTC_adult'
  bulk_samples$cancerType[bulk_samples$cancerType == 'PTC' & bulk_samples$source == 'snRNAseq_Y24'] = 'sn_PTC'
  
  bulk_samples$cancerType[bulk_samples$cancerType == 'PTC' & bulk_samples$source == 'Sanger'] = 'PTC_paediatric'
  bulk_samples$cancerType[bulk_samples$cancerType == 'normal' & bulk_samples$source == 'Sanger'] = 'Normal'
  
  table(bulk_samples$cancerType)
  table(bulk_samples$cancerType,bulk_samples$source)
  
  
  
  
  
  
  
  
  
  
  ##----- 3. Filter out lowly expressed genes -----##
  rawCnt = raw_count
  bulk_dge = DGEList(counts = rawCnt, genes = rownames(rawCnt))
  
  ## Plot library size
  libSize = colSums(rawCnt)
  bulk_samples$libSize = libSize[match(bulk_samples$sampleID,names(libSize))]
  bulk_samples$tech = ifelse(grepl('FFPE',bulk_samples$sampleName),'FFPE','None')
  bulk_samples$source2 = ifelse(grepl('aPTC_|aThy_',bulk_samples$source),'scaThy',bulk_samples$source)
  bulk_samples$cancerType[bulk_samples$source2 == 'scaThy'] = bulk_samples$source[bulk_samples$source2 == 'scaThy']
  p = ggplot(bulk_samples,aes(cancerType,libSize))+
    geom_boxplot(aes(fill=source),outlier.size = 0.001)+
    geom_jitter(aes(color=tech),size=0.4,width = 0.1)+
    facet_grid(.~source2,scales = 'free_x',space = 'free_x')+
    scale_color_manual(values = c('red','black'))+
    theme_classic()+
    #ggtitle(title)+
    scale_y_log10()+
    xlab('')+
    theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
          strip.background=element_rect(linewidth=0),
          axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))
  print(p)
  
  rawCnt_geneDetection = (rawCnt>0)
  rawCnt_geneDetection = colSums(rawCnt_geneDetection)
  bulk_samples$nGene = rawCnt_geneDetection[match(bulk_samples$sampleID,names(rawCnt_geneDetection))]
  p = ggplot(bulk_samples,aes(libSize,nGene,color=tech))+
    geom_point(size=0.3)+
    #geom_boxplot(aes(fill=source),outlier.size = 0.001)+
    #geom_jitter(aes(),size=0.7,width = 0.1)+
    facet_grid(.~source,scales = 'free_x')+
    scale_color_manual(values = c('red','black'))+
    theme_classic()+
    #ggtitle(title)+
    scale_y_log10()+
    theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
          strip.background=element_rect(linewidth=0),
          axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))
  
  print(p)
  
  
  ## Filter out lowly expressed genes
  prop_expressed = rowMeans(edgeR::cpm(bulk_dge,log=TRUE) > 1)
  keep = prop_expressed > 0.1
  op = par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  hist(edgeR::cpm(bulk_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
  abline(v = log(1), lty = 2, col = 2)
  hist(edgeR::cpm(bulk_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
  abline(v = log(1), lty = 2, col = 2)
  par(op)
  
  
  
  ## Subset the count matrix
  bulk_dge = bulk_dge[keep, , keep.lib.sizes = FALSE]
  
  rawCnt = rawCnt[keep, ]
  bulk_dge = DGEList(counts = rawCnt, genes = rownames(rawCnt))
  tpmCnt = tpmCnt[rownames(tpmCnt) %in% rownames(rawCnt),]
  
  cpmCnt = edgeR::cpm(bulk_dge, log = TRUE, prior.count = 1)
  # lapply(split(colnames(cpmCnt),bulk_samples$source),function(e){
  #   dataset = unique(bulk_samples$source[bulk_samples$sampleID %in% e])
  #   print(hist(cpmCnt[,e],main = dataset))})
  
  return(list('rawCnt'=rawCnt,
              'tpmCnt'=tpmCnt,
              'cpmCnt'=cpmCnt,
              'bulk_dge'=bulk_dge,
              'bulk_samples'=bulk_samples
  ))
  
}


##---------------------------------------------##
##   0. Process inhouse bulk RNA-seq data    ####
##---------------------------------------------##
bulk_raw_counts = read.delim('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3426/results/expression_tables/salmon_reads.counts.tsv',sep = '\t')
rownames(bulk_raw_counts) = bulk_raw_counts$ensemblID
bulk_raw_TPM = read.delim('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3426/results/expression_tables/salmon_reads.TPM.tsv',sep = '\t')
rownames(bulk_raw_TPM) = bulk_raw_TPM$ensemblID

row_data = bulk_raw_counts[,c("ensemblID", "geneSymbol", "geneLength","effLength")]

col_data = readxl::read_excel('/lustre/scratch126/cellgen/behjati/mt22/projectManifest.xlsx',sheet = 'bulk_thyroid')
dim(col_data)
col_data$sampleID = paste0('X',col_data$sangerSampleID)
col_data = col_data[col_data$sampleID %in% colnames(bulk_raw_counts),]
dim(col_data)
table(colnames(bulk_raw_counts)[!colnames(bulk_raw_counts) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")] %in% paste0('X',col_data$sangerSampleID))

col_data = as.data.frame(col_data)
rownames(col_data) = col_data$sampleID
col_data$source = 'Sanger'
col_data$cancerType = dplyr::case_when(grepl('Normal',col_data$Tissue) ~ 'Normal',
                                       grepl('Tumour',col_data$Tissue) ~ 'PTC',TRUE ~ 'others')
colnames(col_data)[colnames(col_data) == 'Sex'] = 'sex'
colnames(col_data)[colnames(col_data) == 'Age'] = 'age'
col_data[,c('sampleID','source','sampleName','cancerType','age','sex')]


# Save this as a SummarizedExperiment object
bulk_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts_raw = bulk_raw_counts[,!colnames(bulk_raw_counts) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")],
                counts_tpm = bulk_raw_TPM[,!colnames(bulk_raw_TPM) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")]),
  rowData = row_data,
  colData = col_data,
  metadata = list(genomeVersion = 'GRCh38_v32_2020A-full',
                  gtf_filepath = '/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_modified.gtf')
)
saveRDS(bulk_se,'/lustre/scratch126/cellgen/behjati/mt22/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS')



##-------------------------------------------##
## Module score for fTFC1/fTFC2 signatures ####
##-------------------------------------------##

geneModule = read.csv('/lustre/scratch126/cellgen/behjati/mt22/SupplementaryTableS8_fTFC1.2_geneSignatures.csv')
#colnames(geneModule) = c('Ensembl_ID','Gene_symbol','Chromosome','log2FC fTFC2-vs-fTFC1','logCPM','F','PValue','FDR','percentage cells expressed - fTFC1','percentage cells expressed - fTFC2','DE direction','Gene signature')


## Define gene module
geneList = list('fTFC1' = geneModule$Ensembl_ID[geneModule$DE.direction == 'fTFC2_down'],
                'fTFC2' = geneModule$Ensembl_ID[geneModule$DE.direction == 'fTFC2_up'])



##--------------------------------##
##      Score in bulk data      ####
##--------------------------------##


##--- Define the gene module (convert geneSymbol --> ensID)
moduleList = geneList

moduleList[['fTFC2_combined']] = list('up' = moduleList[['fTFC2']],
                                      'down' = moduleList[['fTFC1']])

moduleList[['fTFC1_combined']] = list('down' = moduleList[['fTFC2']],
                                      'up' = moduleList[['fTFC1']])

##--- import bulk counts and calculate cpmCnt in xx01_moduleScoring.R
bulkRNA = import_bulkRNA_thyroid(bulk_sources = c('inhouse'))
bulk_samples = bulkRNA[['bulk_samples']]
cpmCnt = bulkRNA[['cpmCnt']]
tpmCnt = bulkRNA[['tpmCnt']]
rawCnt = bulkRNA[['rawCnt']]

##---  Score the modules using TPM counts -----##

mtx = tpmCnt[,!colnames(tpmCnt) %in% c('ensID','geneLength')]
# apply the rankGenes method
bulk_ranked = rankGenes(mtx)


# apply the scoring function
allScore = data.frame()
for(i in 1:length(moduleList)){
  if(length(moduleList[[i]]) == 2){
    moduleScores = simpleScore(bulk_ranked,
                               upSet = moduleList[[i]][['up']],
                               downSet = moduleList[[i]][['down']])
    moduleScores = moduleScores[,c('TotalScore', 'TotalDispersion')]
  }else{
    moduleScores = simpleScore(bulk_ranked,upSet = moduleList[[i]])
  }
  
  # create a dataframe with the data required: scores and sample group
  scoredf = merge(bulk_samples,moduleScores,by.x=0,by.y=0)
  scoredf$moduleType = paste0(names(moduleList)[i])
  
  ## Add to allScore
  allScore = rbind(allScore,scoredf)
}


table(allScore$moduleType)
write.csv(allScore,'~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/SingScore_bulkRNA_fTFC1.2_signature_2412.csv')


id_col = 'Sample'
title = 'foetalSig module'

fig4c_fThy_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  
  # create a dataframe with the data required: scores and sample group
  allScore$cancerType[grepl('FFPE',allScore$sampleName)] = paste0('FFPE_',allScore$cancerType[grepl('FFPE',allScore$sampleName)])
  allScore$cancerType[allScore$source == 'Sanger' & allScore$sampleName %in% c('PR66788b','PR66789b','PR66790b','PR66791b','PR66792b')] = 'Normal_foetal'
  allScore$sampleCol = ifelse(grepl('FFPE',allScore$sampleName),'FFPE','normal')
  allScore$sampleCol[allScore$source2 == 'scaThy'] = gsub(':.*$','',allScore$sampleID[allScore$source2 == 'scaThy'])
  allScore$source[allScore$source == 'Yoo_2021'] = 'Yoo_2016'
  allScore$source = factor(allScore$source,c('aPTC_Pu21','aPTC_Wang22','aThy_Hong23','aThy_Mosteiro23',
                                             'GTEx_Thyroid','TCGA_Thyroid','Yoo_2016','He_2021',
                                             'stJudes_Thyroid','Lee_2021','Sanger','scRNAseq_fThy','snRNAseq_Y24.Y46'))
  
  allScore$group_facet_hor = allScore$source
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$cancerType
  
  allScore$group_facet_ver = allScore$moduleType
  #dd = allScore[allScore$age == 'foetus' & allScore$source == 'Sanger' & allScore$moduleType %in% c('fTFC1','fTFC2'),]
  
  # plotFun_sc.fThy.moduleScore_in_Sanger.Fetal.BulkSamples = function(noFrame=FALSE,noPlot=FALSE){
  #   
  #   p1 = ggplot(dd, aes(moduleType, TotalScore)) +
  #     geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.4,linewidth=0.3,fill=grey(0.7)) +
  #     geom_quasirandom(size=0.4,width = 0.15,alpha=0.6)+
  #     scale_y_continuous(breaks = c(0,0.1,0.2,0.3),labels = c(0.0,0.1,0.2,0.3),limits = c(0,0.2))+
  #     theme_classic()+
  #     #ggtitle(title)+
  #     xlab('')+ylab('Module score')+
  #     theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
  #           strip.background=element_rect(linewidth=0),
  #           axis.text = element_text(colour = 'black'),
  #           axis.ticks = element_line(colour = 'black'),
  #           axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))
  #   
  #   print(p1)
  # }
  # 
  # saveFig(file.path(plotDir,'Fig4b_fTFC1.2_moduleScore_bulk.Foetal.Samples'),plotFun_sc.fThy.moduleScore_in_Sanger.Fetal.BulkSamples,rawData=dd,width = 1.6,height = 4,res = 500,useDingbats = F)
  # 
  
  
  plotFun_fTFC2_combined_moduleScore = function(noFrame=FALSE,noPlot=FALSE){
    
    allScore$group_facet_ver = allScore$moduleType
    
    dd = allScore[grepl('FFPE|Thyrocytes|Tumour|fTFC|fThy|Metastatic|Normal|aTFC|C\\d|PTC',allScore$cancerType) & 
                    !grepl('follicular|tallCell|Metastatic|Primary',allScore$cancerType) &
                    #!grepl('Primary',allScore$cancerType_details) &
                    allScore$source %in% c('Sanger','TCGA_Thyroid',#'Yoo_2016',
                                           'He_2021') &
                    allScore$moduleType == 'fTFC2_combined',]
    
    dd$ageGroup = ifelse(dd$source == 'Sanger',dd$ageCat,'adult')
    dd = dd[dd$ageGroup != 'foetus',]
    dd$cancerNormal = ifelse(dd$cancerType %in% c('Normal'),'Normal',
                             ifelse(dd$cancerType == 'Normal.adj','Normal.adj','Tumour'))
    dd$cancerNormal = factor(dd$cancerNormal,c('Normal','Normal.adj','Tumour'))
    dd$med_normal = NA
    for(dataset in unique(dd$source)){
      tmp = dd[dd$source == dataset,]
      med_normal = median(tmp$TotalScore[tmp$cancerNormal == 'Normal'])
      dd$med_normal[dd$source == dataset] = med_normal
    }
    
    dd$normalised_score = dd$TotalScore - dd$med_normal
    dd$source = factor(dd$source,c('Sanger','TCGA_Thyroid','He_2021'))
    
    table(dd$cancerNormal,dd$cancerNormal,dd$source)
    
    
    p1 = ggplot(dd, aes(cancerNormal, normalised_score)) +
      geom_boxplot(aes(fill=cancerNormal),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3,colour='black') +
      geom_quasirandom(size=0.4,width = 0.15,alpha=0.6)+
      scale_fill_manual(values =c(grey(0.8),grey(0.4),'#511378'))+
      #scale_fill_manual(values =c(col25,pal34H))+
      #scale_color_manual(values =c(col25,pal34H))+
      geom_hline(yintercept = 0)+
      #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~source,scales = 'free',space = 'free_x')+
      theme_classic()+
      ylim(-0.1,0.1)+
      #ggtitle(title)+
      xlab('')+ylab('Centralised fTFC2 signature score')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'))
    
    print(p1)
  }
  # saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 7,height = 5,res = 500,useDingbats = T)
  # saveFig(file.path(plotDir,'Fig3c_fTFC1.2_moduleScore_bulkSamples_sub'),plotFun_fTFC2_combined_moduleScore,rawData=allScore,width = 4.8,height = 4,res = 500,useDingbats = F)
  
}