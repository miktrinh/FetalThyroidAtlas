library(Seurat)  #updates from 3.2.3 to 4.0.3
library(tidyverse)
source('~/lustre_mt22/generalScripts/utils/misc.R')
# 
# library(patchwork)
# library(cowplot)
# library(ggplot2)
# library(Signac)
# library(tidyverse)
# library(org.Hs.eg.db)
# library(clusterProfiler)
# library(psych)
# library(qgraph)
# library(igraph)
# library(GSVA)
# library(GSEABase)
# library(limma)
# library(hdf5r)


#input h5 data
T1L.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743021_T1L.h5", use.names = T)
T1R.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743022_T1R.h5", use.names = T)
T2L.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743023_T2L.h5", use.names = T)
T2R.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743024_T2R.h5", use.names = T)
T3L.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743025_T3L.h5", use.names = T)
T3R.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743026_T3R.h5", use.names = T)
NT.data <- Read10X_h5("~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/GSM5743027_NT.h5", use.names = T)

# Create Seurat object and quality control
T1L <- CreateSeuratObject(counts = T1L.data, project = "Thyroid_L1", min.cells = 3, min.features = 200)
T1L$stim <- "LEFT1"
T1L[["percent.mt"]] <- PercentageFeatureSet(T1L, pattern = "^MT-")
#VlnPlot(T1L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T1L <- subset(T1L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T2L <- CreateSeuratObject(counts = T2L.data, project = "Thyroid_L2", min.cells = 3, min.features = 200)
T2L$stim <- "LEFT2"
T2L[["percent.mt"]] <- PercentageFeatureSet(T2L, pattern = "^MT-")
#VlnPlot(T2L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2L <- subset(T2L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T3L <- CreateSeuratObject(counts = T3L.data, project = "Thyroid_L3", min.cells = 3, min.features = 200)
T3L$stim <- "LEFT3"
T3L[["percent.mt"]] <- PercentageFeatureSet(T3L, pattern = "^MT-")
#VlnPlot(T3L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T3L <- subset(T3L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

NT <- CreateSeuratObject(counts = NT.data, project = "Thyroid_L4", min.cells = 3, min.features = 200)
NT$stim <- "Normal"
NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")
#VlnPlot(NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NT <- subset(NT, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T1R <- CreateSeuratObject(counts = T1R.data, project = "Thyroid_R1", min.cells = 3, min.features = 200)
T1R$stim <- "RIGHT1"
T1R[["percent.mt"]] <- PercentageFeatureSet(T1R, pattern = "^MT-")
#VlnPlot(T1R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T1R <- subset(T1R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T2R <- CreateSeuratObject(counts = T2R.data, project = "Thyroid_R2", min.cells = 3, min.features = 200)
T2R$stim <- "RIGHT2"
T2R[["percent.mt"]] <- PercentageFeatureSet(T2R, pattern = "^MT-")
#VlnPlot(T2R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2R <- subset(T2R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

T3R <- CreateSeuratObject(counts = T3R.data, project = "Thyroid_R3", min.cells = 3, min.features = 200)
T3R$stim <- "T3R"
T3R[["percent.mt"]] <- PercentageFeatureSet(T3R, pattern = "^MT-")
#VlnPlot(T3R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T3R <- subset(T3R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)


##Perform integration using CCA method

thyroid.anchors <- FindIntegrationAnchors(object.list = list(T1L, T2L, T3L, NT, T1R, T2R, T3R), dims = 1:20)
thyroid.combined <- IntegrateData(anchorset = thyroid.anchors, dims = 1:20)
thyroid.combined

#Perform an integrated analysis

DefaultAssay(thyroid.combined) <- "integrated"
#NormalizeData
thyroid.combined <- NormalizeData(thyroid.combined, normalization.method = "LogNormalize", scale.factor = 10000)
#FindVariableFeatures
thyroid.combined <- FindVariableFeatures(thyroid.combined, selection.method = "vst", nfeatures = 2000)

#to check nFeature_RNA", "nCount_RNA", "percent.mt
VlnPlot(thyroid.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Run the standard workflow for visualization and clustering
thyroid.combined <- ScaleData(thyroid.combined, verbose = FALSE)
thyroid.combined <- RunPCA(thyroid.combined, npcs = 30, verbose = FALSE)


#Determine the 'dimensionality' of the dataset 
thyroid.combined <- JackStraw(thyroid.combined, num.replicate = 100)
thyroid.combined <- ScoreJackStraw(thyroid.combined, dims = 1:20)
JackStrawPlot(thyroid.combined, dims = 1:15)
ElbowPlot(thyroid.combined)


# t-SNE and Clustering
thyroid.combined <- RunUMAP(thyroid.combined, reduction = "pca", dims = 1:15)
thyroid.combined <- FindNeighbors(thyroid.combined, reduction = "pca", dims = 1:15)
thyroid.combined <- FindClusters(thyroid.combined, resolution = 0.8)
thyroid.combined <- RunTSNE(thyroid.combined, reduction = "pca", dims = 1:15)
thyroid.combined$cellID = colnames(thyroid.combined)
thyroid.combined$sampleID = thyroid.combined$stim
thyroid.combined@meta.data$stim <- factor(thyroid.combined@meta.data$stim, 
                                          levels = c("LEFT1", "LEFT2", "LEFT3", "Normal", "RIGHT1", "RIGHT2", "T3R"), 
                                          labels = c("T1L", "T2L", "T3L", "NT", "T1R", "T2R", "T3R"))

DimPlot(thyroid.combined, reduction = "umap", group.by = "sampleID",cols = col25)
DimPlot(thyroid.combined, reduction = "umap", group.by = "seurat_clusters",cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()

#name the clusters
thyroid.combined <- RenameIdents(thyroid.combined, `0` = "Follicular_cells", `2` = "Follicular_cells",`4` = "Follicular_cells", `5` = "Follicular_cells",`9` = "Follicular_cells", `7` = "Follicular_cells",`11` = "Follicular_cells",`22` = "Follicular_cells",`24` = "Follicular_cells",  
                                 `1` = "T_cells", `3` = "T_cells", `10` = "T_cells",`13` = "T_cells",`17` = "T_cells",`15` = "T_cells",
                                 `8` = "Myeloid_cells", `27` = "Myeloid_cells", 
                                 `14` = "Endothelial_cells", `26` = "Endothelial_cells", `25` = "Endothelial_cells", `31` = "Endothelial_cells", `23` = "Endothelial_cells", `30` = "Endothelial_cells", `16` = "Endothelial_cells", 
                                 `18` = "Pericyte",`21` = "Pericyte",`6` = "Pericyte",`12` = "Pericyte",`20` = "Fibroblast",
                                 `19` = "B_cells",`28` = "B_cells", `29` = "B_cells")
thyroid.combined$celltype = as.character(Idents(thyroid.combined))

View(thyroid.combined@meta.data)


# Visualization for Figure 1C and 1D
DimPlot(thyroid.combined, reduction = "umap",group.by = 'celltype',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()


#VlnPlot for Figure 1E
DefaultAssay(thyroid.combined) = 'RNA'
thyroid.combined <- NormalizeData(thyroid.combined)
thyroid.combined <- ScaleData(thyroid.combined)
VlnPlot(thyroid.combined, features = c("TG", "HIGD1B","CD3D", "CD68", "VWF", "CD79A","PDGFRA","TPSAB1"), stack = T, flip = T)


# find markers for every cluster compared to all remaining cells, report only the positive ones
thyroid.combined.markers <- FindAllMarkers(thyroid.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(thyroid.combined.markers,"Supplementary_File_1_DEGs_among_8_known_cell_clusters.csv")

#heatmap Figure 1F
DoHeatmap(thyroid.combined, features = c("TG","KRT18","KRT19","TSHR","KRT7"
                                         ,"HIGD1B","CSRP2","CACNB2","COL25A1","RGS5"
                                         ,"CD3G","GNLY","PDCD1","CD8B","FOXP3"
                                         ,"CD86","LYZ","HLA-DRA","APOC1","S100A8"
                                         ,"VWF","ARL15","PLPP1","PTPRG","STC1"
                                         ,"CD79A","DERL3","BANK1","MS4A1","LY9"
                                         ,"PDGFRA","COL3A1","COL1A1","MGP","COL1A2"
                                         ,"KIT","CSF1","S100B"), angle = 90) + NoLegend()


##Fraction of cell types in each sample for Figure 1G
table(thyroid.combined$stim,thyroid.combined$celltype)
write.csv(table(thyroid.combined$stim,thyroid.combined$celltype),"celltypefrequency.csv")

## Add metadata
thyroid.combined$mutation = ifelse(thyroid.combined$stim == 'T1L','RET_FARP1_fusion',
                                   ifelse(thyroid.combined$stim == 'NT','Normal','BRAF_V600E'))
thyroid.combined$PTC_subtype = ifelse(thyroid.combined$stim %in% c('T2L','T2R'),'Follicular',
                                      ifelse(thyroid.combined$stim %in% c('NT'),'Normal','Classical'))

mdat = cbind(thyroid.combined@meta.data,as.data.frame(thyroid.combined@reductions$umap@cell.embeddings))
write.csv(mdat,'~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/Wang_etal_2022_mdat.csv')
saveRDS(thyroid.combined,'~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/Wang_etal_2022.RDS')



##----------------------------##
##   Set Global parameters  ####
##----------------------------##
library(GenomicFeatures)
source('~/lustre_mt22/generalScripts/utils/pseudobulk.R')

#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

## Generic gene map
geneMap = read.table('~/lustre_mt22/Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236839_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',sep = '\t')
colnames(geneMap) = c('ensID','geneSym','library')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns)[match(geneMap$ensID,gns$gene_id)])




##-------------------------##
##  aPTC vs Normal aTFC  ####
##-------------------------##
thyroid.combined = readRDS('~/lustre_mt22/Thyroid/Data/published_scRNAseq/wang_etal_2022/Wang_etal_2022.RDS')
DimPlot(thyroid.combined,group.by = 'mutation',cols = col25,label = T,repel = T,label.box = T)

thyroid.combined$group = ifelse(thyroid.combined$celltype == 'Follicular_cells',thyroid.combined$mutation,thyroid.combined$celltype)
DimPlot(thyroid.combined,cells.highlight = thyroid.combined$cellID[thyroid.combined$group == 'BRAF_V600E'])

Idents(thyroid.combined) = thyroid.combined$group

## Import pPTC markers
tum.vs.norm_markers.sub = read.csv('pThy_tum.vs.norm_markers_topMarkers_240519.csv')
tum.vs.norm_markers.sub = tum.vs.norm_markers.sub[order(abs(tum.vs.norm_markers.sub$pct.diff),decreasing = T),]

ret_markers = FindMarkers(thyroid.combined,ident.1 = 'RET_FARP1_fusion',ident.2='Normal',features = tum.vs.norm_markers.sub$gene)
ret_markers$gene = rownames(ret_markers)
ret_markers$ensID = geneMap$ensID[match(ret_markers$gene,geneMap$geneSym)]
rownames(ret_markers) = ret_markers$ensID
ret_markers = annotateGenes(ret_markers,geneMap = geneMap)

#write.csv(ret_markers,'~/lustre_mt22/Thyroid/Results_v2/04_pThyCancer_Tumour.vs.Thyrocytes/aThy.Wang22_RET.tum.vs.norm_markers_allGenes_2405.csv')
write.csv(ret_markers,'~/lustre_mt22/Thyroid/Results_v2/04_pThyCancer_Tumour.vs.Thyrocytes/aThy.Wang22_RET.tum.vs.norm_markers_pPTC.degs_240527.csv')

ret_markers = ret_markers[ret_markers$p_val_adj < 0.05,]
ret_markers$direction = ifelse(ret_markers$avg_log2FC > 0,'tum_up','tum_down')
ret_markers$pct.diff = ret_markers$pct.1 - ret_markers$pct.2
ret_markers.sub = ret_markers[abs(ret_markers$avg_log2FC) > 1 & abs(ret_markers$pct.diff) > 0.2,]

thyroid.combined$group2 = ifelse(thyroid.combined$celltype == 'Follicular_cells' & thyroid.combined$mutation == 'Normal','Normal',
                                 ifelse(thyroid.combined$celltype == 'Follicular_cells' & thyroid.combined$mutation != 'Normal','Tumour',thyroid.combined$celltype))
Idents(thyroid.combined) = thyroid.combined$group2
aTum_markers = FindMarkers(thyroid.combined,ident.1 = 'Tumour',ident.2='Normal')
aTum_markers$gene = rownames(aTum_markers)
aTum_markers$ensID = geneMap$ensID[match(aTum_markers$gene,geneMap$geneSym)]
rownames(aTum_markers) = aTum_markers$ensID
aTum_markers = annotateGenes(aTum_markers,geneMap = geneMap)
write.csv(aTum_markers,'~/lustre_mt22/Thyroid/Results_v2/04_pThyCancer_Tumour.vs.Thyrocytes/aThy.Wang22_allTum.vs.norm_markers_allGenes_2405.csv')

aTum_markers = aTum_markers[aTum_markers$p_val_adj < 0.05,]
aTum_markers$direction = ifelse(aTum_markers$avg_log2FC > 0,'tum_up','tum_down')
aTum_markers$pct.diff = aTum_markers$pct.1 - aTum_markers$pct.2
aTum_markers.sub = aTum_markers[abs(aTum_markers$avg_log2FC) > 1 & abs(aTum_markers$pct.diff) > 0.2,]




##----    Import pThy tum vs norm markers ------##
pTum_markers = read.csv('~/lustre_mt22/Thyroid/Results_v2/04_pThyCancer_Tumour.vs.Thyrocytes/pThy_tum.vs.norm_markers_allGenes_240514.csv')
pTum_markers = pTum_markers[pTum_markers$p_val_adj < 0.05,]
pTum_markers$direction = ifelse(pTum_markers$avg_log2FC > 0,'tum_up','tum_down')
pTum_markers$pct.diff = pTum_markers$pct.1 - pTum_markers$pct.2
pTum_markers.sub = pTum_markers[abs(pTum_markers$avg_log2FC) > 1 & abs(pTum_markers$pct.diff) > 0.2,]


library(UpSetR)
upset(fromList(list(aTum_up = ret_markers$gene[ret_markers$direction == 'tum_up'],
                    aTum_down = ret_markers$gene[ret_markers$direction == 'tum_down'],
                    pTum_up = pTum_markers$geneSym[pTum_markers$direction == 'tum_up'],
                    pTum_down = pTum_markers$geneSym[pTum_markers$direction == 'tum_down'])))


upset(fromList(list(aTum_up = ret_markers.sub$gene[ret_markers.sub$direction == 'tum_up'],
                    aTum_down = ret_markers.sub$gene[ret_markers.sub$direction == 'tum_down'],
                    pTum_up = pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up'],
                    pTum_down = pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_down'])))

upset(fromList(list(aTum_up = aTum_markers.sub$gene[aTum_markers.sub$direction == 'tum_up'],
                    aTum_down = aTum_markers.sub$gene[aTum_markers.sub$direction == 'tum_down'],
                    pTum_up = pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up'],
                    pTum_down = pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_down'])))



genes_toPlot = pTum_markers.sub$gene[pTum_markers.sub$direction == 'tum_up' ]
genes_toPlot = pTum_markers.sub$gene[pTum_markers.sub$direction == 'tum_up' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_up']]
genes_toPlot = ret_markers.sub$gene[ret_markers.sub$direction == 'tum_up' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up']]

genes_toPlot = pTum_markers.sub$gene[pTum_markers.sub$direction == 'tum_down' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_down']]
genes_toPlot = ret_markers.sub$gene[ret_markers.sub$direction == 'tum_down' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_down']]


##---- in aPTC
Idents(thyroid.combined) = thyroid.combined$group
DotPlot(thyroid.combined,
        #features = 
        features = genes_toPlot
        ) +
  
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('')


##---- in pPTC
pPTC = readRDS('~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS')
pPTC$group = pPTC$finalAnn_broad
pPTC$group[pPTC$group == 'Tumour'] = ifelse(pPTC$etiology[pPTC$group == 'Tumour'] == 'Met', 'Met:Y46', paste0(pPTC$group[pPTC$group == 'Tumour'],':',pPTC$donor[pPTC$group == 'Tumour']))
pPTC$group[pPTC$group == 'Thyrocytes'] = ifelse(pPTC$etiology[pPTC$group == 'Thyrocytes'] == 'Met', 'Met_Normal:Y46', paste0(pPTC$group[pPTC$group == 'Thyrocytes'],':',pPTC$donor[pPTC$group == 'Thyrocytes']))
pPTC$group = factor(pPTC$group,c('Thyrocytes:Y24','Thyrocytes:Y46','Met_Normal:Y46','Tumour:Y24','Tumour:Y46','Met:Y46',
                                 unique(pPTC$group[!grepl('Thyrocytes|Normal|Tumour|Met',pPTC$group)])))

Idents(pPTC) = pPTC$group
DotPlot(pPTC,idents = c('Thyrocytes:Y24','Thyrocytes:Y46','Met_Normal:Y46','Tumour:Y24','Tumour:Y46','Met:Y46'),
        #features = ret_markers.sub$gene[!is.na(ret_markers.sub$ensID) & ret_markers.sub$direction == 'tum_up' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up']]
        features = genes_toPlot
        
        ) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('')

##---- in combined fThy/aThy
combinedSrat = readRDS('~/lustre_mt22/Thyroid/Results_v2/x04_adult.vs.foetal_thyrocytes/fThy_aThy1_aThy2_combinedSrat.RDS')
combinedSrat$finalAnn_lvl2 = ifelse(combinedSrat$finalAnn %in% c('aTFC3','thy_Lumen-forming'),'TFC2',
                                    ifelse(combinedSrat$finalAnn %in% c('aTFC1','aTFC2','aTFC4','aTFC5','C0','C1','C2','C3','thy_TH_processing'),'TFC1',
                                           ifelse(combinedSrat$finalAnn %in% c('C4'),'C4','others')))
combinedSrat$pcw[is.na(combinedSrat$pcw)] = 'adult'

combinedSrat$group = paste0(combinedSrat$pcw,':',combinedSrat$finalAnn)
Idents(combinedSrat) = combinedSrat$group
#DotPlot(combinedSrat,idents = unique(combinedSrat$finalAnn[grepl('TFC|thy_|C\\d',combinedSrat$finalAnn)]),

DotPlot(combinedSrat,idents = unique(combinedSrat$group[grepl('thy_Lu|thy_TH',combinedSrat$group)]),group.by = 'pcw',
        features = genes_toPlot
)+RotatedAxis()+
  theme(axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')




library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
selected_attributes = c('ensembl_gene_id','external_gene_name','chromosome_name','description','gene_biotype')
aThy_Wang22_geneMap = getBM(attributes = selected_attributes,values = rownames(thyroid.combined),filters = 'external_gene_name',mart = ensembl)



moduleList = list('unique_pPTC_up' = pTum_markers.sub$ensID[pTum_markers.sub$direction == 'tum_up' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_up']],
                  'unique_aPTC_up' = ret_markers.sub$ensID[!is.na(ret_markers.sub$ensID) & ret_markers.sub$direction == 'tum_up' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up']],
                  'unique_pPTC_down' = pTum_markers.sub$ensID[pTum_markers.sub$direction == 'tum_down' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_down']],
                  'unique_aPTC_down' = ret_markers.sub$ensID[!is.na(ret_markers.sub$ensID) & ret_markers.sub$direction == 'tum_down' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_down']])


moduleList = list('unique_pPTC' = list('up'=pTum_markers.sub$ensID[pTum_markers.sub$direction == 'tum_up' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_up']],
                                       'down'=pTum_markers.sub$ensID[pTum_markers.sub$direction == 'tum_down' & !pTum_markers.sub$gene %in% ret_markers.sub$gene[ret_markers.sub$direction == 'tum_down']]),
                  'unique_aPTC' = list('up'=ret_markers.sub$ensID[!is.na(ret_markers.sub$ensID) & ret_markers.sub$direction == 'tum_up' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_up']],
                                       'down' = ret_markers.sub$ensID[!is.na(ret_markers.sub$ensID) & ret_markers.sub$direction == 'tum_down' & !ret_markers.sub$gene %in% pTum_markers.sub$geneSym[pTum_markers.sub$direction == 'tum_down']]))
                  
