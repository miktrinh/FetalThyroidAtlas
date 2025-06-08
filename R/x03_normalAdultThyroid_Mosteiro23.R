## process Normal adult thyroid from Mosteiro et.al. 2023
## then run LR using REF = foetal thyroid

library(tidyverse)
library(Seurat)
library(Matrix)
source('/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_utils.R')
source('/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R')

# dataDir = '/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/aThyroid_Mosteiro_2023/aThy_Mosteiro2023_matrix/'
# if(!dir.exists(dataDir)){
#   print('Cannot find dataDir, please check!')
# }
# 
# srat = Read10X(dataDir)
# 
# colnames(srat)
# mdat = read.csv(file.path(dataDir,'metadata.csv'))
# colnames(mdat)[colnames(mdat) == 'X'] = 'cellID'
# rownames(mdat) = mdat$cellID
# mdat = mdat[match(colnames(srat),mdat$cellID),]
# 
# srat = CreateSeuratObject(srat,meta.data = mdat)
# srat = standard_clustering(srat)
# 
# umap = read.csv(file.path(dataDir,'UMAP.csv'),row.names = 1)
# umap = as.matrix(umap[match(colnames(srat),mdat$cellID),])
# rownames(umap) = mdat$cellID
# colnames(umap) = c('UMAP_1','UMAP_2')
# srat@reductions$umap@cell.embeddings = umap
# 
# ## Add metadata
# srat$annot = ifelse(srat$leiden_scvi == 0,'aTFC1',
#                     ifelse(srat$leiden_scvi == 1,'aTFC2',
#                            ifelse(srat$leiden_scvi == 2,'aTFC3',
#                                   ifelse(srat$leiden_scvi == 4,'aTFC4',
#                                          ifelse(srat$leiden_scvi == 5,'aTFC5',
#                                                 ifelse(srat$leiden_scvi == 3,'aEndo/Mes','others'))))))
# View(srat@meta.data)
# DimPlot(srat,group.by = 'annot',label = T,label.box = T,repel = T,cols = col25[-6]) + NoLegend()
# 
# saveRDS(srat,'/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/aThyroid_Mosteiro_2023/aThyroid_Mosteiro_2023_sratObj.RDS')
# 
# srat = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/aThyroid_Mosteiro_2023/aThyroid_Mosteiro_2023_sratObj.RDS')
# 
# library(SoupX)
# qm = quickMarkers(srat@assays$RNA@counts,srat$annot)


# #-------------------------------------------------------------------------------------------##
# 
# ## LR similarity between fThy Thyrocytes and paediatric thyrocytes
# 
# #------------------#
# ##    Library     ##
# #------------------#
# library(tidyverse)
# library(Seurat)
# library(RColorBrewer)
# library(readxl)
# 
# ##---------     HASSAN'S LR METHOD - L2 (similar to celltypist)   ---------####
# ## NOTE: cannot reproduce Hassan's result - need to ask him how to get 0.95 probability ....
# result = read.csv('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/x_LR_hm11/fThy2n.REF_aMosteiro23.tgt_predictions.csv')
# result$annot = srat$annot[match(result$index,srat$cellID)]
# 
# ggplot(result,aes(annot,probabilities,fill=celltype_predictions))+
#   geom_boxplot(outlier.size = 0.0001)+
#   theme_classic()+theme(legend.position = 'none')
# View(result[result$annot == 'aTFC1',])
# 
# ## Heatmap of mean_probability per cluster
# meanProb = do.call(rbind,lapply(split(result,result$annot),function(i){
#   annot = unique(i$annot)
#   df = i[,grepl('^predict_',colnames(i))]
#   df = t(colMeans(df)) %>% as.data.frame()
#   df$annot = annot
#   return(df)
# }))
# 
# 
# medianProb = do.call(rbind,lapply(split(result,result$annot),function(i){
#   annot = unique(i$annot)
#   df = i[,grepl('^predict_',colnames(i))]
#   df = apply(df,2,function(x){return(median(x,na.rm = T))}) %>% t() %>%  as.data.frame()
#   df$annot = annot
#   return(df)
# }))
# 
# library(ComplexHeatmap)
# mtx = meanProb[,!grepl('annot',colnames(meanProb))] %>% as.matrix()
# mtx = meanProb[,grepl('thy_',colnames(meanProb))] %>% as.matrix()
# Heatmap(((t(mtx))),show_row_names = T,show_column_names = T)








##-----------------------------------------------------------------------------------------##
####   3. LRv2_wCT with REF = foetal Thyrocytes (from Hassan) ; tgt = aThy_Mosteiro_23   ####
##-----------------------------------------------------------------------------------------##
## AIM: Train LRv2_wCT model on foetal thyrocytes (REF)
##      Predict similarity on: 1. aThy 2. aThy - Mosteiro 3. pPTC


source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/sharedCode-main/logisticRegressionCellTypist.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_utils.R")


outDir = '/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/5_LRv2_fThy2nREF/mar24'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)











##----------------------------------------##
##    1. Train model - REF = fThy2n     ####
##----------------------------------------##
## Import fThy2n object
REF.srat = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/fetalThyroid/fThyroid_annotated_fromHassan_jul23.RDS')
REF.srat = FindVariableFeatures(REF.srat,nfeatures = 5000)
## Train LRv2 model
message(sprintf('Training new model on srat object of %d cells.',ncol(REF.srat)))

model_fp = '/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/5_LRv2_fThy2nREF/mar24/fThy2nREF_hvg_LRv2wCT_model.pkl'

if(!file.exists(model_fp)){
  trainCelltypistModel(cnts = REF.srat@assays$RNA@counts[rownames(REF.srat@assays$RNA@counts) %in% REF.srat@assays$RNA@var.features,],
                       outPath = model_fp,
                       labels = REF.srat$celltype,
                       n_jobs = -1)
  
}else{
  message('Existing model found. Skipping...')
}

message('1. completed model training on REF.srat')


##------------------------------##
##    2. Prepare tgt.srat     ####
##        aThy_v1               ##
##        aThy_v2               ##
##        pPTC                  ##
##------------------------------##
tgt.srat_fp = file.path(outDir,'tgtSrat_aThy1_aThy2_pPTC.RDS')
if(file.exists(tgt.srat_fp)){
  tgt.srat = readRDS(tgt.srat_fp)
}else{
  ## Import adult srat object
  aThy_1 = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/Hong_etal_2023/GSE182416_sratObj.RDS')
  
  ## Import adult Mosteiro 2023 object
  aThy_2 = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/aThyroid_Mosteiro_2023/aThyroid_Mosteiro_2023_sratObj.RDS')
  
  ## Import pPTC
  pPTC = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/1_thyroid_annotation/thyroid_clean_ann_soupedRho0.2_jul23.RDS')
  pPTC$annot_lvl1 = as.character(pPTC$finalAnn)
  pPTC$annot_lvl2 = as.character(pPTC$finalAnn)
  pPTC$annot_lvl2[grepl('^thy',pPTC$annot_lvl1)] = 'Thyrocytes'
  pPTC$annot_lvl2[grepl('pDC|MonoMac|basophil',pPTC$annot_lvl1)] = 'Monocytes'
  pPTC$annot_lvl2[grepl('NK|T\\.cell',pPTC$annot_lvl1)] = 'T_cells'
  pPTC$annot_lvl2[grepl('\\?',pPTC$annot_lvl1)] = 'Fibroblast'
  pPTC$annot_lvl2[pPTC$seurat_clusters == 57] = '57'
  pPTC$annot_lvl2[pPTC$annot_lvl2 == 'Thyrocytes' & pPTC$etiology == 'left_inferior_tumour'] = 'Tumour'
  #DimPlot(pPTC,group.by = 'annot_lvl2',label = T,label.box = T,repel = T) + NoLegend() 
  
  
  ## Combine all 3 objects into tgt.srat
  tgt.srat = merge_seurat_objects(aThy_1,aThy_2,keepAllGenes = F,genomeVersions = c('v38','v38'))
  tgt.srat = merge_seurat_objects(tgt.srat,pPTC,keepAllGenes = F,genomeVersions = c('v38','v38'))
  #tgt.srat = standard_clustering(tgt.srat)
  tgt.srat$dataset = ifelse(tgt.srat$cellID %in% aThy_2$cellID,'aThy_Mosteiro',
                            ifelse(tgt.srat$cellID %in% pPTC$cellID,'pPTC','aThy_1'))
  tgt.srat$finalAnn = ifelse(tgt.srat$dataset == 'aThy_Mosteiro',tgt.srat$annot,
                             ifelse(tgt.srat$dataset == 'pPTC',tgt.srat$annot_lvl2,tgt.srat$celltype_global))
  tgt.srat$finalAnn[tgt.srat$finalAnn == 'Thyrocytes' & tgt.srat$etiology == 'left_inferior_tumour'] = 'Tumour'
  
  saveRDS(tgt.srat,file.path(outDir,'tgtSrat_aThy1_aThy2_pPTC.RDS'))
  
}

#DimPlot(tgt.srat,group.by = 'dataset',cols=c(col25,pal34H))

##-----------------------------------------##
##    3. Predict similarity on tgt.srat  ####
##-----------------------------------------##
ref_tissue = 'thyroid'
message(sprintf('Predicting label against fetal %s',ref_tissue))

skipIfExists = T
lr_output_fp = file.path(outDir,'fThy.REF_onto_aThy.pPTC_LRwCT_output.RDS')
if(file.exists(lr_output_fp) & skipIfExists){
  lr_output = readRDS(lr_output_fp)
}else{
  if(!dir.exists(dirname(lr_output_fp))){
    dir.create(dirname(lr_output_fp))
  }
  
  lr_output = runCelltypist(cnts=tgt.srat@assays$RNA@counts,
                            model=model_fp)
  saveRDS(lr_output,lr_output_fp)
  
  message('Prediction completed!')
  
}  






# # Remove rubbish genes
# geneExpr_pThy = rowSums(pThy@assays$RNA@counts)
# geneExpr_aThy = rowSums(aThy@assays$RNA@counts)
# geneExpr_fThy = rowSums(REF.srat@assays$RNA@counts)
# 
# genesToKeep = geneExpr_pThy[geneExpr_pThy>=10 & names(geneExpr_pThy) %in% names(geneExpr_aThy[geneExpr_aThy >= 10]) & 
#                               names(geneExpr_pThy) %in% names(geneExpr_fThy[geneExpr_fThy >= 10])]
# 
# genesToKeep = genesToKeep[names(genesToKeep) %in% rownames(REF.srat)]
# genesToKeep = genesToKeep[names(genesToKeep) %in% rownames(aThy)]
# genesToKeep = names(genesToKeep)
# genesToKeep = genesToKeep[!grepl('^MT-|^RPL|RPS|MALAT1|NEAT1|AC\\d+|AL\\d+',genesToKeep)]
# thy_prog = c('PAX8','NKX2-1','FOXE1','HHEX','UROD','DIO2','DUOX1','DUOX2','DUOXA1','DUOXA2',
#              'SLC5A5','ANO1','SLC26A4','TPO','IYD','TG','PAX8','GLIS3','TSHR','SLC16A2','SLC16A10','EXOC4','ELMO1','VPS13C','STON2','SPG11')
# # Remove soupy genes
# soupyGenes = read.csv('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/fThy2nREF_downsampled_Thyrocytes_markersExprinpThy.csv')
# soupyGenes = soupyGenes$gene[soupyGenes$frac_cellExpr >= 0.5]
# genesToKeep = genesToKeep[!genesToKeep %in% c(soupyGenes,thy_prog,'PAX8','TSHR','EEF1A1','SFTA3')]
# 
# 
# 
# pThy@assays$RNA@counts = pThy@assays$RNA@counts[genesToKeep,]
# aThy@assays$RNA@counts = aThy@assays$RNA@counts[genesToKeep,]
# REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
# 
# ## Merge pThy and aThy
# tgt.srat = merge_seurat_objects(srat1 = pThy, srat2 = aThy,keepAllGenes = F,genomeVersions = c('v38','v38'))
# tgt.srat$type = ifelse(tgt.srat$cellID %in% aThy$cellID,'aThyroid',
#                        ifelse(tgt.srat$cellID %in% pThy$cellID,'pThyroid','others'))
# 
# REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[rownames(tgt.srat),]
# 
# print(table(rownames(REF.srat) %in% rownames(tgt.srat)))
# message(('2. tgt.srat loaded'))
# 
# 
# 
# 
# ## Configure REF and tgt sratObj for LR
# # Only keep the same genes between REF and tgt.srat
# REF.srat$type = 'fThyrocytes'
# tgt.srat$donor[is.na(tgt.srat$donor)] = tgt.srat$orig.ident_new[is.na(tgt.srat$donor)]
# tgt.srat$type = ifelse(tgt.srat$donor == 'Y24', 'pThyroid','aThyroid')
# tgt.srat$annot = ifelse(is.na(tgt.srat$finalAnn),tgt.srat$celltype_sub,tgt.srat$finalAnn)
# tgt.srat$seurat_clusters = paste0(tgt.srat$type,':',as.character(tgt.srat$seurat_clusters))
# 
# 
# 
# #####################################################################################
# ##----------------------------------------##
# ####    3. Run LR_orig on tgt.srat      ####
# ##----------------------------------------##
# 
# 
# ##----------------------------##
# ##   Run Logistic Regression  ##
# ##----------------------------##
# skipIfExists=F
# ref_annot='annot'
# maxCells=4000
# #outDir = file.path('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/1_thyroid_annotation/LRorig_fThyrocyteREF_jul23/')
# model_fp = '/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/adult_and_paediatric/highLevel_REFlabel/nov23/fThy2nREF_trainModel_4kmaxcells_70perc_geneFiltered_highLevel.RDS'
# 
# LR_level='both'
# srat_annot=''
# minGeneMatch = 0.99
# maxCells=4000
# tissue = 'thyroid'
# 
# #out_prefix = 'fThyroidREF.highLevel_on_p.aThyrocyte_maxCells_70perc_'
# #plot_prefix = 'fThyroidREF.highLevel_on_p.aThyrocyte_maxCells_70perc_'
# out_prefix = 'fThy2nREF.highLevel_on_pThy.aThy_maxCells_70perc_'
# plot_prefix = 'fThy2nREF.highLevel_on_pThy.aThy_maxCells_70perc_'
# 
# plot_prefix = NULL
# 
# outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
#                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
#                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
#                 scLR_TGTtype='',scLR_REFtype='')
# 
# message(sprintf('3. LR completed for tissue %s',tissue))
# 
# 
# message(sprintf('Saving temp. annotated tgt.srat@meta.data object for aThyroid'))
# write.csv(tgt.srat@meta.data, file.path(outDir,'p.aThyroid.highLevelLR_mdat_231128.csv'))
# 
# ## Import output
# outputs = readRDS(file.path(outDir,'fThy2nREF.highLevel_on_pThy.aThy_maxCells_70perc_raw_LR_outputs.RDS'))
# 
# 
# 
# 
# 
# ##---------------------------------##
# ##   Do some LR_similarity plots   ##
# ##---------------------------------##
# ## annotated Cluster level #
# # if(length(outputs) > 2){
# #   output = outputs[[1]]
# # }else{
# #   output = outputs[[2]][[1]]
# # }
# # 
# # type = ifelse(grepl('ref_',rownames(output)),'REF',
# #               ifelse(rownames(output) %in% tgt.srat$seurat_clusters[tgt.srat$type == 'pThyroid'],'pThyroid','aThyroid'))
# # show_row_names = T
# # column_order = colnames(output)[order(colnames(output))]
# # row_order = rownames(output)[!grepl('Tumour',rownames(output))]
# # row_order = row_order[order(row_order)]
# # row_order = c(row_order,rownames(output)[grepl('Tumour',rownames(output))])
# # 
# # 
# # plot_prefix = 'fThyroidref.highLevel_on_p.aThyroid_maxCells_70perc_'
# # 
# # pdf(file.path(outDir,paste0(plot_prefix,'clusterLR.pdf')),width = 12,height = 15)
# # 
# # hm = similarityHeatmap(output,
# #                        row_order=row_order,
# #                        column_order = column_order,
# #                        row_title_rot = 0,
# #                        row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
# # draw(hm)
# # dev.off()
# # 
# # 
# # ## annotated Single-cell level #
# # if(length(outputs) > 2){
# #   output = outputs[[2]]
# # }else{
# #   output = outputs[[2]][[2]]
# # }
# # 
# # 
# # #in_mtx = output[rownames(output) %in% c(tgt.srat$cellID[tgt.srat$type == 'pThyroid']),]
# # in_mtx = output[rownames(output) %in% c(tgt.srat$cellID[tgt.srat$annot %in% c('C0','C1','C2','C3','C4','Cytotoxic CD8T','Fibroblasts',
# #                                                                               'imm_T.cell','mes_Fibroblast','Plasma','pDC','smc_Pericyte','thy_Thyrocytes_LRRK2high','thy_Thyrocytes_LRRK2low','thy_Thyrocytes_mt')],
# #                                         paste0('ref_',REF.srat$cellID[grepl('^thy_|smc_Pericytes|imm_T_cells',REF.srat$annot)])),]
# # output = in_mtx
# # tgt.srat$annot2 = ifelse(tgt.srat$type == 'pThyroid' & tgt.srat$etiology == 'left_inferior_tumour',
# #                          paste0('p:T:',tgt.srat$annot),
# #                          ifelse(tgt.srat$type == 'pThyroid' & tgt.srat$etiology != 'left_inferior_tumour',
# #                                 paste0('p:N:',tgt.srat$annot),paste0('a:',tgt.srat$annot)))
# # type = ifelse(grepl('ref_',rownames(output)),
# #               as.character(REF.srat$annot[match(rownames(output),paste0('ref_',REF.srat$cellID))]),
# #               as.character(tgt.srat$annot2[match(rownames(output),tgt.srat$cellID)]))
# # 
# # 
# # 
# # df = data.frame(cellID = rownames(output),annot = type)
# # m = match(df$cellID,tgt.srat$cellID)
# # df$annot.orig = tgt.srat$annot2[m]
# # table(df$annot == df$annot.orig)
# # 
# # type[grepl('^p:T:thy_',type)] = 'p:T:Thyrocytes'
# # type[grepl('^p:N:thy_',type)] = 'p:N:Thyrocytes'
# # type = factor(type,levels = c(unique(type[grepl('^a',type)]),
# #                               unique(type[grepl('^p',type) & !grepl('^p:.+:Thyrocytes',type)]),
# #                               unique(type[grepl('^p:.+:Thyrocytes',type)]),
# #                               unique(type[!grepl('^a|^p',type)])))
# # pdf(file.path(outDir,paste0(plot_prefix,'scLR_REF.subset.pdf')),width = 8,height = 20)
# # show_row_names=F
# # hm = similarityHeatmap(output,
# #                        column_order = colnames(output)[order(colnames(output))],
# #                        row_title_rot = 0,
# #                        row_title_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=12),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=12),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
# # draw(hm)
# # 
# # dev.off()
# # 
# # ## annotated Single-cell level #
# # # output = outputs[[2]][['scLR_all']]
# # # in_mtx = output
# # # type = as.character(tgt.srat$seurat_clusters[match(rownames(output),tgt.srat$cellID)])
# # # type[is.na(type)] = as.character(REF.srat$annot[match(rownames(output)[is.na(type)],paste0('ref_',REF.srat$cellID))])
# # # type = factor(type, levels = c(colnames(output)[order(colnames(output))],seq(0:max(as.numeric(tgt.srat$seurat_clusters)))))
# # # pdf(file.path(outDir,paste0(plot_prefix,'scLR_all.pdf')),width = 6,height = 30)
# # # show_row_names=F
# # # hm = similarityHeatmap(in_mtx,
# # #                        column_order = colnames(output)[order(colnames(output))],
# # #                        row_title_rot = 0,
# # #                        row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
# # #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# # #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
# # # draw(hm)
# # # 
# # # dev.off()
# # 
# # 
# # Extract markers of the LR model
# # model = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/adult_and_paediatric/highLevel_REFlabel/nov23/fThy2nREF_trainModel_4kmaxcells_70perc_geneFiltered_highLevel.RDS')
# # markers = getMarkers(model)
# # # Subset for thyrocytes markers
# # markers = markers[grepl('thy_',markers$class),]
# # 
# # # plot some top markers for each fThy population
# # markers = markers[order(abs(markers$coef),decreasing = T),]
# # 
# # ## Extract only thyrocytes from reference foetal thyroid object
# # DotPlot(thy,group.by = 'annot',scale = T,
# #         features = unique(c(#markers$gene[markers$coef > 0.25 & markers$class == 'thy_Cycling'],
# #           markers$gene[markers$coef > 0.05 & markers$class == 'thy_Lumen-forming'],
# #           markers$gene[markers$coef > 0.01 & markers$class == 'thy_TH_processing']
# #           #markers$gene[markers$coef > 0.01 & markers$class == 'imm_Monocytes']
# #         )))+RotatedAxis()+theme(axis.text.x = element_text(size = 8))
# # 
# # 
# # DotPlot(pThy,group.by = 'etiology',scale = T,
# #         features = unique(c(#markers$gene[markers$coef > 0.25 & markers$class == 'thy_Cycling'],
# #           markers$gene[markers$coef > 0.05 & markers$class == 'thy_Lumen-forming'],
# #           markers$gene[markers$coef > 0.01 & markers$class == 'thy_TH_processing']
# #         )))+RotatedAxis()+theme(axis.text.x = element_text(size = 8))
# # 
# # DotPlot(aThy,group.by = 'celltype_sub',scale = T,
# #         features = unique(c(#markers$gene[markers$coef > 0.25 & markers$class == 'thy_Cycling'],
# #           markers$gene[markers$coef > 0.05 & markers$class == 'thy_Lumen-forming'],
# #           markers$gene[markers$coef > 0.01 & markers$class == 'thy_TH_processing']
# #         )))+RotatedAxis()+theme(axis.text.x = element_text(size = 8))
# # 
# # 
# # DotPlot(aThy,group.by = 'celltype_sub',scale = F,
# #         features = unique(c(markers$gene[markers$coef > 0.25 & markers$class == 'thy_Cycling'],
# #                             markers$gene[markers$coef > 0.05 & markers$class == 'thy_Lumen-forming'],
# #                             markers$gene[markers$coef > 0.01 & markers$class == 'thy_TH_processing']
# #         )))+RotatedAxis()
# # FeaturePlot(aThy,'TPO')
# # 
# # 
# # harm.f.a.Thyroid = merge_seurat_objects(srat1 = REF.srat,srat2 = aThy,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # harm.f.a.Thyroid$annot = ifelse(!is.na(harm.f.a.Thyroid$celltype),harm.f.a.Thyroid$celltype,harm.f.a.Thyroid$celltype_sub)
# # harm.f.a.Thyroid$batch = ifelse(is.na(harm.f.a.Thyroid$donor),harm.f.a.Thyroid$orig.ident_new,harm.f.a.Thyroid$donor)
# # harm.f.a.Thyroid = standard_clustering(harm.f.a.Thyroid,runHarmony = T,harmonyVar = c('batch'))
# # harm.f.a.Thyroid$dataset = ifelse(grepl('^Old|^Young',harm.f.a.Thyroid$batch),'aThyroid','fThyroid')
# # DimPlot(harm.f.a.Thyroid,group.by = 'annot',label = T,repel = T,label.box = T,cols = col25) +NoLegend()
# # DimPlot(harm.f.a.Thyroid,group.by = 'dataset',label = T,repel = T,label.box = T,cols = col25)
# # 
# # 
# # 
# # 
# # 
# # ## EnrichR on thy_Procesing markers
# # ## Enrichment analysis ###
# # library(enrichR)
# # listEnrichrSites()
# # setEnrichrSite("Enrichr") # Human genes
# # websiteLive <- TRUE
# # dbs <- listEnrichrDbs()
# # if (is.null(dbs)) websiteLive <- FALSE
# # if (websiteLive) head(dbs)
# # 
# # dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021",
# #          'HDSigDB_Human_2021','KEGG_2021_Human','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
# #          'TF_Perturbations_Followed_by_Expression','WikiPathways_2019_Human',
# #          'Cancer_Cell_Line_Encyclopedia','CCLE_Proteomics_2020','CellMarker_Augmented_2021',
# #          'Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up',
# #          'Drug_Perturbations_from_GEO_2014','Drug_Perturbations_from_GEO_down','Drug_Perturbations_from_GEO_up','DrugMatrix','IDG_Drug_Targets_2022',
# #          'Elsevier_Pathway_Collection','Enrichr_Submissions_TF-Gene_Coocurrence','Gene_Perturbations_from_GEO_down','Gene_Perturbations_from_GEO_up')
# # 
# # thyProc_markers_enrichR <- enrichr(markers$gene[markers$class == 'thy_TH_processing' & markers$gene > 0.01], dbs)
# # thyLumen_markers_enrichR <- enrichr(markers$gene[markers$class == 'thy_Lumen-forming' & markers$gene > 0.01], dbs)
# # thyCycling_markers_enrichR <- enrichr(markers$gene[markers$class == 'thy_Cycling' & markers$gene > 0.01], dbs)
# # plotEnrich(thyProc_markers_enrichR[[5]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# # 
# # View(thyProc_markers_enrichR[[5]])
# # 
# # 
# # 
# # 
