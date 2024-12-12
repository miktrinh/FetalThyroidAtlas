

library(readxl)
library(SummarizedExperiment)
library(edgeR)
library(singscore)
library(GenomicFeatures)
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

import_bulkRNA_thyroid = function(bulk_sources = c('Lee2021','Yoo2016','TCGA_Thyroid','snPaedThyroid','scaThy','GTEx_Thyroid','StJudes_Thyroid','fAdr','He2021','inhouse')){
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
  
  
  if('He2021' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/He_etal_2021/aPTC_He_2021_sce.RDS'
    he21_sce = readRDS(sce_path)
    sce_list[['He2021']] = he21_sce
    he21_mdat = as.data.frame(colData(he21_sce))
    
    bulk_samples = rbind(bulk_samples,he21_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    he21_rawCnt = assays(he21_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = he21_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(he21_rawCnt))
      raw_count = cbind(raw_count[genesToKeep,],he21_rawCnt[genesToKeep,])
    }
  }
  
  if('Lee2021' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/Lee_etal_2021/pedPTC_Lee_2021_sce.RDS'
    pedPTC_sce = readRDS(sce_path)
    sce_list[['Lee2021']] = pedPTC_sce
    pedPTC_mdat = as.data.frame(colData(pedPTC_sce))
    
    pedPTC_mdat$sampleName = pedPTC_mdat$Sample.Name
    pedPTC_mdat$age = pedPTC_mdat$Age
    bulk_samples = rbind(bulk_samples,pedPTC_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    pedPTC_rawCnt = assays(pedPTC_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = pedPTC_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(pedPTC_rawCnt))
      raw_count = cbind(raw_count[genesToKeep,],pedPTC_rawCnt[genesToKeep,])
    }
    
  }
  
  
  if('Yoo2016' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/Yoo_etal_2016/adultPTC_Yoo_2016_sce.RDS'
    aThy_sce = readRDS(sce_path)
    sce_list[['Yoo2016']] = aThy_sce
    
    aThy_mdat = as.data.frame(colData(aThy_sce))
    bulk_samples = rbind(bulk_samples,aThy_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    aThy_rawCnt = assays(aThy_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = aThy_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(aThy_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],aThy_rawCnt[genesToKeep,])
    }
    
  }
  
  if('TCGA_Thyroid' %in% bulk_sources){
    tcga_sce_path_m2 = '~/lustre_mt22/Thyroid/Data/TCGA_Thyroid/TCGA_Thyroid_gdc0923_sce.RDS'
    tcga_sce = readRDS(tcga_sce_path_m2)
    sce_list[['TCGA_Thyroid']] = tcga_sce
    tcga_mdat = as.data.frame(colData(tcga_sce))
    
    bulk_samples = rbind(bulk_samples,tcga_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    tcga_rawCnt = assays(tcga_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = tcga_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(tcga_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],tcga_rawCnt[genesToKeep,])
    }
    
  }
  
  
  if('snPaedThyroid' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/paedThyroid_snRNAseq_pseudobulk_sce.RDS'
    pThyroid_sce = readRDS(sce_path)
    sce_list[['snPaedThyroid']] = pThyroid_sce
    
    pThyroid_mdat = as.data.frame(colData(pThyroid_sce))
    bulk_samples = rbind(bulk_samples,pThyroid_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    pThyroid_rawCnt = assays(pThyroid_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = pThyroid_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(pThyroid_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],pThyroid_rawCnt[genesToKeep,])
    }
  }
  
  if('scfThy' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/fThy_scRNAseq_pseudobulk_sce.RDS'
    scfThy_sce = readRDS(sce_path)
    sce_list[['scfThy']] = scfThy_sce
    
    scfThy_mdat = as.data.frame(colData(scfThy_sce))
    bulk_samples = rbind(bulk_samples,scfThy_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    scfThy_rawCnt = assays(scfThy_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = scfThy_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(scfThy_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],scfThy_rawCnt[genesToKeep,])
    }
  }
  
  
  
  if('scaThy' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/aThy_scRNAseq_pseudobulk_sce.RDS'
    scaThy_sce = readRDS(sce_path)
    sce_list[['scaThy']] = scaThy_sce
    
    scaThy_mdat = as.data.frame(colData(scaThy_sce))
    bulk_samples = rbind(bulk_samples,scaThy_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    scaThy_rawCnt = assays(scaThy_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = scaThy_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(scaThy_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],scaThy_rawCnt[genesToKeep,])
    }
  }
  
  
  if('GTEx_Thyroid' %in% bulk_sources){
    gtex_sce_path = '~/lustre_mt22/Thyroid/Data/GTEx_Thyroid/GTEx_Thyroid_sce.RDS'
    gtex_sce = readRDS(gtex_sce_path)
    sce_list[['GTEx_Thyroid']] = gtex_sce
    
    gtex_mdat = as.data.frame(colData(gtex_sce))
    bulk_samples = rbind(bulk_samples,gtex_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    gtex_rawCnt = assays(gtex_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = gtex_rawCnt
    }else{
      rowDat = rowData(gtex_sce)
      #rownames(gtex_rawCnt) = rowDat$ensID[match(rownames(gtex_rawCnt),rowDat$gene_id)]
      genesToKeep = intersect(rownames(raw_count),rownames(gtex_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],gtex_rawCnt[genesToKeep,])
    }
  }
  
  
  if('StJudes_Thyroid' %in% bulk_sources){
    stJudes_sce_path = '~/lustre_mt22/Thyroid/Data/StJudes_Thyroid/StJudes_Thyroid_230921_sce.RDS'
    stJudes_sce = readRDS(stJudes_sce_path)
    sce_list[['StJudes_Thyroid']] = stJudes_sce
    
    stJudes_mdat = as.data.frame(colData(stJudes_sce))
    bulk_samples = rbind(bulk_samples,stJudes_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    stJudes_rawCnt = assays(stJudes_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = stJudes_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(stJudes_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],stJudes_rawCnt[genesToKeep,])
    }
  }
  
  if('inhouse' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/inhouse_bulkRNA_thyroid/inhouse_bulkRNA_thyroid_2410_sce.RDS'
    inhouse_sce = readRDS(sce_path)
    sce_list[['Sanger']] = inhouse_sce
    inhouse_mdat = as.data.frame(colData(inhouse_sce))
    
    bulk_samples = rbind(bulk_samples,inhouse_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    inhouse_rawCnt = assays(inhouse_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = inhouse_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(inhouse_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],inhouse_rawCnt[genesToKeep,])
    }
  }
  
  if('fAdr' %in% bulk_sources){
    sce_path = '~/lustre_mt22/Thyroid/Data/fAdr.SCPs_scRNAseq_pseudobulk_sce.RDS'
    scfAdr_sce = readRDS(sce_path)
    sce_list[['fAdr']] = scfAdr_sce
    
    scfAdr_mdat = as.data.frame(colData(scfAdr_sce))
    bulk_samples = rbind(bulk_samples,scfAdr_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    
    scfAdr_rawCnt = assays(scfAdr_sce)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = scfAdr_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(scfAdr_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],scfAdr_rawCnt[genesToKeep,])
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





ggRLE <- function(dat_x, annot, col_str,isLog=TRUE, isLarge=FALSE,colorVal=col25,
                  ylim = c(-2,2),zero_line=TRUE, zero_col="skyblue", medPoint=FALSE, whisk=1.5){
  
  col_title <- str_to_title(col_str)
  
  if (!(setequal(colnames(dat_x), rownames(annot)))){
    message("Make sure annot rownames match dat_x sample names")
  }
  
  if (!(isLog)){
    dat_x <- log2(dat_x+1)
  }
  
  ## Subset to common sample names
  annot <- annot[order(as.vector(annot[[col_str]])), , drop = F]
  annot$ColourBy <- as.vector(annot[[col_str]])
  annot$Sample <- rownames(annot)
  dat_x <- dat_x[, rownames(annot)]
  
  ## RLE boxplots
  rle <- dat_x - rowMedians(dat_x)
  rleLong <- reshape2::melt(rle, value.name = "RLE", varnames = c("genes", "Sample"))
  rleLong$Sample <- as.character(rleLong$Sample)
  
  rleLong <- merge(data.table::data.table(varhandle::unfactor(rleLong)),  # faster merging
                   data.table::data.table(varhandle::unfactor(annot)),
                   by = "Sample", sort=F)
  
  ## Calculate the median of the RLE boxplots:
  rleLong <- rleLong %>%
    group_by(Sample) %>%
    mutate(MedRLE = median(RLE)) %>%
    ungroup() %>%
    data.frame()
  rleLong$Sample <- factor(rleLong$Sample , levels=unique(rleLong$Sample))
  
  ## Remove Whiskers of the boxplots if the sample size is very large -- Sep's code
  if(isLarge){
    whisk=0
  }
  
  if (medPoint){
    gg <- ggplot(rleLong, aes(x = Sample, y = RLE, fill = ColourBy))+
      stat_boxplot(geom = "errorbar", width = 0.3)+
      geom_boxplot(outlier.shape = NA, coef=whisk)+
      geom_point(data = rleLong[! duplicated(rleLong$Sample), ],
                 aes(x = Sample, y = MedRLE, fill = ColourBy),
                 size = 2, shape = 21, colour = "black", lwd = 2)+
      scale_fill_manual(name=col_title,values = colorVal)+
      scale_y_continuous(name = "RLE",limits = ylim)+theme_bw()+
      theme(panel.grid.minor = element_blank(),
            #  panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    #axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=8,colour = 'black'))
  } else {
    gg <- ggplot(rleLong, aes(x = Sample, y = RLE, fill = ColourBy))+
      stat_boxplot(geom = "errorbar", width = 0.3)+
      geom_boxplot(outlier.shape = NA,coef=whisk)+
      scale_fill_manual(name=col_title,values = colorVal)+
      scale_y_continuous(name = "RLE",limits = ylim)+theme_bw()+
      theme(panel.grid.minor = element_blank(),
            #  panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    #axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=8,colour = 'black'))
  }
  
  if (zero_line){
    gg <- gg + geom_hline(yintercept = 0, col = zero_col, lwd = 1)+
      theme(panel.grid.major=element_blank())
  }
  return(gg)
}

