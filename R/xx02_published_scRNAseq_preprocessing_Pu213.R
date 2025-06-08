## Process aPTC from Pu etal 2021

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")


dataDir = list.dirs('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/GSE184362_RAW',full.names = T,recursive = F)
names(dataDir) = basename(dataDir)
names(dataDir) = gsub('_','.',names(dataDir))
dataDir = dataDir[!names(dataDir) %in% c('.ipynb.checkpoints','cleanCounts')]

srat = Read10X(dataDir)
srat = CreateSeuratObject(srat)
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

srat = subset(srat,subset = nCount_RNA >= 500 & nFeature_RNA >= 300 & percent.mt <= 30)
srat = standard_clustering(srat)
DimPlot(srat,group.by = 'orig.ident',cols = col25,label = T,repel = T,label.box = T)

FeaturePlot(srat,c('TG','TPO'))

##--------------------------------------##
##    LR_v1; REF = fThy; tgt = Pu21   ####
##--------------------------------------##
source("~/lustre_mt22/generalScripts/utils/logisticRegression.R")
# Using same model as LR for PTC
model_fp = '~/lustre_mt22/Thyroid/Results_v2/05_LRv1_fThy.REF/LRv1_fThy.full.REF_pPTC.aPTC.tgt_alpha.0.5_maxCell.4k_maxGeneFilter/LRv1_fThy.full.REF_trainModel_alpha.0.5_maxCell.4k_maxGeneFilter.RDS'
model = readRDS(model_fp)

# Predict similarity
logit_mtx = predictSimilarity(model, srat@assays$RNA@counts,logits = T,minGeneMatch = 0.95)
# pp = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = F)

# Do heatmap
in_mtx = logit_mtx

type = srat$seurat_clusters[match(rownames(in_mtx),rownames(srat@meta.data))]
table(is.na(type))

pdf('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/LRv1_fThy.full.REF_alpha.0.5_maxCell.4k_maxGeneFilter_scLR.pdf',width = 10,height = 50)
show_row_names=F
hm = similarityHeatmap(in_mtx,
                       column_order = colnames(in_mtx)[order(colnames(in_mtx))],
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
ht = draw(hm)

dev.off()


## Add annotation based on LR
annot = do.call(rbind,apply(in_mtx,1,function(x){
  i = which(x==max(x))
  tmp = data.frame(max_LR = max(x),celltype = colnames(in_mtx)[i])
  return(tmp)
  }))
srat$cellID = rownames(srat@meta.data)
srat@meta.data = cbind(srat@meta.data,annot[match(srat$cellID,rownames(annot)),])
srat$celltype[srat$max_LR < 1] = 'unknown'
DimPlot(srat,group.by = 'celltype',cols = col25)
        #cols = sample(c(col25,pal34H),n_distinct(srat$celltype)))
DimPlot(srat,cells.highlight = srat$cellID[srat$celltype == 'fTFC1' & srat$max_LR < 2])
DimPlot(srat,group.by = 'seurat_clusters',label = T)
saveRDS(srat,'~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/Pu21_annotated_sratObj.RDS')

## Sub clustering thyrocytes/Tumour only
srat.thy = subset(srat, subset = celltype %in% c('fTFC1','fTFC2'))
srat.thy = standard_clustering(srat.thy)
DimPlot(srat.thy,group.by = 'orig.ident',cols = col25)
DimPlot(srat.thy,cells.highlight = srat.thy$cellID[srat.thy$orig.ident == 'PTC9.P'])
DimPlot(srat.thy,cells.highlight = srat$cellID[srat$seurat_clusters == 2 & srat$celltype=='fTFC1'])
DimPlot(srat.thy,cells.highlight = srat$cellID[srat$seurat_clusters == 2 & srat$celltype=='fTFC2'])
srat.thy$group = srat$group[match(srat.thy$cellID,srat$cellID)]
Idents(srat.thy) = srat.thy$group

thyroid_markers = c('TSHR','PAX8','GLIS3','TG', # follicular cells markers
                    'NKX2-1','IYD','TPO','HHEX','FOXE1','DUOXA1','DUOXA2','DUOX1','DUOX2','DIO',
                    'SLC26A4','ANO','SLC5A5','SLC16A2','SLC16A10',
                    'DPP6','ZNF804B','KCNQ1' # thyroid epithelial cells
)

DotPlot(srat.thy,scale = T,idents = c('fTFC1:2','fTFC2:2'),
        features = thyroid_markers
        #features = genes_toPlot
)+RotatedAxis()+
  theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=11),
        #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = 'bottom') + xlab('') + ylab('')

srat$group = ifelse(srat$celltype %in% c('fTFC1','fTFC2') & srat$seurat_clusters == 2, paste0(srat$celltype,':2'),as.character(srat$celltype))
srat$group = as.factor(srat$group)
Idents(srat) = srat$group
ncell_perCT = table(srat$group)
DotPlot(srat,scale = T,idents = names(ncell_perCT[ncell_perCT > 30 & !names(ncell_perCT) %in% c('fTFC1','fTFC2')]),
        features = thyroid_markers
        #features = genes_toPlot
)+RotatedAxis()+
  theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=11),
        #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = 'bottom') + xlab('') + ylab('')






##----------------------------------------------------------------##
##    LR_v1; REF = fTFC1 / fTFC2 / SCPs; tgt = Pu21_normalThy   ####
##----------------------------------------------------------------##
source("~/lustre_mt22/generalScripts/utils/logisticRegression.R")
# Using same model as LR for PTC
model_fp = '~/lustre_mt22/Thyroid/Results_v2/05_LRv1_fThy.REF/LRv1_fThy.REF_pPTC.aPTC.tgt_alpha.0.5_maxCell.4k_maxGeneFilter/LRv1_fThy.REF_trainModel_alpha.0.5_maxCell.4k_maxGeneFilter.RDS'
model = readRDS(model_fp)

# Predict similarity
logit_mtx = predictSimilarity(model, srat@assays$RNA@counts,logits = T)
# pp = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = F)

# Do heatmap
in_mtx = logit_mtx

type = srat$celltype[match(rownames(in_mtx),rownames(srat@meta.data))]
table(is.na(type))

pdf(file.path(outDir,paste0(plot_prefix,'scLR_all.pdf')),width = 10,height = 50)
pdf('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/LRv1_fThy.REF_alpha.0.5_maxCell.4k_maxGeneFilter_scLR.pdf',width = 10,height = 50)
show_row_names=F
hm = similarityHeatmap(in_mtx,
                       column_order = colnames(in_mtx)[order(colnames(in_mtx))],
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
ht = draw(hm)

dev.off()

##--- Ellie's style plot
use_logit = F
mtx = in_mtx
if(!use_logit){
  mtx = (1+exp(-mtx))**-1
}

mtx = as.data.frame(mtx)
mtx$cellID = rownames(mtx)
df = pivot_longer(mtx,cols = colnames(mtx)[colnames(mtx) != 'cellID'],names_to = 'REF_celltype',values_to = 'LRpp_score')
df$group = srat$celltype[match(df$cellID,srat$cellID)]

if(use_logit){
  # df$LR_score_group = ifelse(df$LRpp_score > 3, 'high',
  #                            ifelse(df$LRpp_score > 0, 'mid2',
  #                                   ifelse(df$LRpp_score < -3,'low','mid')))
  
}else{
  df$LR_score_group = ifelse(df$LRpp_score >= 0.9, 'high',
                             ifelse(df$LRpp_score >= 0.6, 'mid2',
                                    ifelse(df$LRpp_score >= 0.4, 'white',
                                           ifelse(df$LRpp_score <= 0.1,'low','mid'))))  
}


df$LR_score_group = factor(df$LR_score_group,levels = rev(c('high','mid2','white','mid','low')))

df$toKeep = ifelse(grepl('fTFC',df$group),T,F)
data_toPlot = df[df$toKeep == T,]
data_toPlot$REF_celltype = factor(data_toPlot$REF_celltype,levels = c('fTFC1','fTFC2','SCPs','B_cells'))

# source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
# plotDir='~/lustre_mt22/Thyroid/Figures'

plotFun_LRv1 = function(noFrame=FALSE,noPlot=FALSE){
  p1=ggplot(data_toPlot,aes(y=group,fill = LR_score_group))+
    geom_bar(position = 'fill',col='black',linewidth=0.3,width=0.7)+
    facet_grid(.~(REF_celltype),scales = 'free_y',space = 'free')+
    #scale_fill_manual(values = rev(c('#EA2335','#F06571','white','#2D3F90'))) +
    scale_fill_manual(values = rev(c('#EA2335','#F7B2B8','white','#B3B9DC','#2D3F90'))) +
    scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1))+
    theme_classic(base_size = 13) + xlab('Fraction of cells') + ylab('')+
    theme(#panel.border = element_rect(fill=F,color='black',linewidth = 1),
      panel.border = element_blank(),
      axis.line.y = element_blank(),axis.ticks.y = element_blank(),
      strip.background = element_blank(),axis.text = element_text(colour = 'black'))
  
  print(p1)
}


DimPlot(srat.thy,cells.highlight = srat$cellID[srat$celltype == 'fTFC2' & srat$seurat_clusters == 2])


