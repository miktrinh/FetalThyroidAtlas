library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

# Converting between Mouse-Human gene list
mouse_hum_bioMartconverter = function(){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mouse_genelist=getBM(mart = mouse,attributes = c('mgi_symbol'))
  homGenes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.vector(mouse_genelist$mgi_symbol), 
                    mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  write_delim(genesV2, '/nfs/team292/mt22/homGemes_mmus_hsap.txt', delim = '\t',col_names = T)
  
}

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x, 
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(humanx)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , 
                   mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(humanx)
}


inhouse_seurat =function(mtx, meta_data, mt_pattern="^MT-"){
  srat=CreateSeuratObject(mtx,min.cells = 3, min.features  = 200, meta.data = meta_data)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = mt_pattern)
  srat=subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)
  #srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  #srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  #exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
  #keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  #srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  #srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
  
}


standard_clustering = function(srat,nPCs=75,clusteringRes=1,skipCCS=FALSE,s.genes=NULL,g2m.genes=NULL,runHarmony=F,harmonyVar = NULL,doPlot=T,...){
  require(Seurat)
  require(cowplot)
  srat = NormalizeData(srat)
  srat = FindVariableFeatures(srat)
  srat = ScaleData(srat,...)
  srat = RunPCA(srat, npcs = nPCs)
  
  if(runHarmony){
    # check the harmony variables provided
    if(is.null(harmonyVar) | !all(harmonyVar %in% colnames(srat@meta.data))){
      stop('No or not all variables for Harmony provided exist...')
    }else{
      require(harmony)
      
      srat = RunHarmony(srat,harmonyVar, plot_convergence = TRUE)
      # To directly access the new Harmony embeddings, use the Embeddings command.
      harmony_embeddings <- Embeddings(srat, 'harmony')
    
      if(doPlot){
        options(repr.plot.height = 5, repr.plot.width = 12)
        p1 <- DimPlot(object = srat, reduction = "harmony", pt.size = .1, group.by = harmonyVar[1])
        p2 <- VlnPlot(object = srat, features = "harmony_1", group.by = harmonyVar[1], pt.size = .1)
        plot_grid(p1,p2)  
      }
    }
  }
  #ElbowPlot(srat, ndims = nPCs)

  # Always cluster without Harmony first
  srat = FindNeighbors(srat, dims=seq(nPCs))
  srat = FindClusters(srat, resolution = clusteringRes)
  srat = RunUMAP(srat, dims=seq(nPCs))
  
  if(!runHarmony | is.null(harmonyVar)){doPlot=F}
  if(doPlot){
    p = list()
    for(var in harmonyVar){
      p[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Pre-Harmony')
    }
    plot_grid(plotlist = p)
  }
  
  if(runHarmony){
    srat = FindNeighbors(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = FindClusters(srat, resolution = clusteringRes, reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = RunUMAP(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    
    if(doPlot){
      pHarm = list()
      for(var in harmonyVar){
        pHarm[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Post-Harmony')
      }
      plot_grid(plotlist = pHarm)
    }
  }
  


  # CellCycle Scoring
  if(!skipCCS){
    if(is.null(s.genes)){
      s.genes <- cc.genes$s.genes  
    }
    if(is.null(g2m.genes)){
      g2m.genes <- cc.genes$g2m.genes  
    }
    
    srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)  
  }
  
  return(srat)
}


pseudobulk_normalization = function(cnt_mtx,genes_to_include,split_by){
  pb = do.call(cbind,lapply(split(colnames(cnt_mtx),split_by),function(e) rowSums(cnt_mtx[,e,drop=FALSE])))
  pb = as.matrix(pb)
  # minibulk_normalization: expression of each gene in each minibulk is normalized (divided by) the total count in the minibulk
  pb_norm = sweep(pb,2,colSums(pb),`/`)
  # log transformation
  pb_norm = log(pb_norm * 10^4 + 1)
  pb_norm = pb_norm[rowSums(pb_norm) != 0,]
  
  # Subset the count matrix to only include cell types of interest
  pb.sub = pb_norm[rownames(pb_norm) %in% genes_to_include,]
  
  # Row normalization (Zscore method by (x - mu_i)/sd_i )
  pb.sub_rowNorm = sweep(pb.sub,1,rowMeans(pb.sub),`-`)
  pb.sub_rowNorm = sweep(pb.sub_rowNorm,1,apply(pb.sub,1,FUN = sd),`/`)
  pb.sub_rowNorm = pb.sub_rowNorm[!is.na(rownames(pb.sub_rowNorm)),]
  
  return(pb.sub_rowNorm)
}



avgExpr_cellFrac_byGroup = function(srat,group=NULL,genes=NULL,doPlot=F){
  if(is.null(genes)){
    genes = rownames(srat)
  }
  if(is.null(group)){
    group = 'seurat_clusters'
  }
  srat$group_tmp = srat@meta.data[[group]]
  out_df = data.frame()
  for(g in unique(srat$group_tmp)){
    if(length(srat$cellID[srat$group_tmp == g]) == 1){next}
    mtx = srat@assays$RNA@counts[genes,srat$cellID[srat$group_tmp == g]]
    tmp = as.data.frame(rowSums(mtx > 0))
    colnames(tmp) = 'nCell'
    tmp$gene = rownames(tmp)
    tmp$totalCell = ncol(mtx)
    tmp$group = g  
    out_df = rbind(out_df,tmp)
  }
  out_df$frac = out_df$nCell / out_df$totalCell
  
  ## Avergae expression
  avgExpr = AverageExpression(srat,features = genes,group.by = 'group_tmp')
  avgExpr_df = as.data.frame(avgExpr$RNA)
  avgExpr_df$gene = rownames(avgExpr_df)
  avgExpr_df = pivot_longer(avgExpr_df,cols = 1:(ncol(avgExpr_df) - 1),names_to = 'group',values_to = 'avgExpr')
  
  table(out_df$group %in% avgExpr_df$group)
  if(doPlot){
    if(sum(out_df$group %in% avgExpr_df$group) > 1){
      out_df_merged = merge(out_df,avgExpr_df,by = c('gene','group'),all=T)
      
      p = ggplot(out_df_merged,aes(gene,frac,fill=avgExpr))+
        geom_col()+
        facet_wrap(vars(group),ncol = 4)+
        scale_fill_gradient(low = colAlpha('#3373C4',0.3),high = '#03254C')+
        ylab('Fraction of cells') + xlab('')+
        theme_classic(base_size = 13)+theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
                                            panel.border = element_rect(fill = F,linewidth=1,colour = 'black'),
                                            axis.line = element_blank(),legend.position = 'top')  
      print(p)
    }
  }
  
  
  return(list(cellFrac_df = out_df,
              avgExpr_df = avgExpr_df))
  
}

seurat_integration = function(object.list, resolution, fvf.nfeatures=2000, all_genes = TRUE){
  
  #### Perform integration ####
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = object.list,fvf.nfeatures = fvf.nfeatures)
  anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
  # this command creates an 'integrated' data assay
  if (all_genes == TRUE){
    all_gene_names = c()
    for(i in 1:length(object.list)){
      all_gene_names = unique(c(all_gene_names,rownames(object.list[[i]])))
    }
  }
  if(length(all_gene_names) == 0){
    print('ERROR: all gene names were not retrieved!')
    return()
  }
  combined <- IntegrateData(anchorset = anchors,features.to.integrate = all_gene_names)
  
  #### Integrated Analysis
  # specify that we will perform downstream analysis on the corrected data note that the original
  # unmodified data still resides in the 'RNA' assay
  DefaultAssay(combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  #combined = ScaleData(object = combined, vars.to.regress = c("nCount_RNA", "percent.mt",'S.Score','G2M.Score'))
  combined = ScaleData(object = combined)
  combined <- RunPCA(combined, npcs = 60, verbose = FALSE)
  #ElbowPlot(combined, ndims = 60)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined,resolution = resolution)
  
  return(combined)
}


## Most of the time, Seurat::merge() is fine. But because it keeps the union of genes between the 2 objects
# this might be problematic if the 2 objects were generated and pre-processed differently (eg. published external data or different genome versions etc.)
# also, should always match on geneID instead of gene names

merge_seurat_objects = function(srat1, srat2, keepAllGenes=F, genomeVersions = NULL){
  require(tidyverse)
  ## Check genome version used to map data from each of the 2 objects
  if(is.null(genomeVersions) | length(genomeVersions) != 2){
    stop('Please provide genome versions used for the 2 seurat objects')
  }
  
  ## Map gene names via ensID
  if(genomeVersions[1] == genomeVersions[2]){
    message(sprintf('Same genome version (%s) was used for both objects - No further geneID processing is needed',unique(genomeVersions)))
  }else if(genomeVersions[1] != genomeVersions[2]){
    message(sprintf('Mapping geneID for both objects to v38'))
  }
  
  ## use Seurat function to generate merged metadata
  srat1@meta.data$cellID = rownames(srat1@meta.data)
  srat2@meta.data$cellID = rownames(srat2@meta.data)
  merged.mdat = plyr::rbind.fill(srat1@meta.data,srat2@meta.data)
  
  if(sum(duplicated(merged.mdat$cellID)) > 0){
    rownames(merged.mdat) =  make.names(merged.mdat$cellID, unique=TRUE)  
  }else{
    rownames(merged.mdat) = merged.mdat$cellID
  }
  
  #rownames(merged.mdat) =  gsub('^X','',rownames(merged.mdat))
  #sum(duplicated(colnames(srat2)))
  
  
  # If we wish to check that this is equivalent to what Seurat::merge() does to merge metadata, here's the code for it
  # I've checked though and its good!
  #require(Seurat)
  #merged.srat = merge(srat1,srat2)
  #merged.mdat2 = merged.srat@meta.data
  #merged.mdat2$cellID = rownames(merged.mdat2)
  #indx <- sapply(merged.mdat, is.factor)
  #merged.mdat[indx] <- lapply(merged.mdat[indx], function(x) as.character(x))
  #dplyr::all_equal(merged.mdat,merged.mdat2)
  require(Seurat)
  
  if(keepAllGenes){
    # Use Seurat::merge() method
    merged.srat = merge(srat1,srat2)
    return(merged.srat)
  }
  
  ## If NOT keepAllGenes - ie. only keep genes that are common between the 2 srat objects. 
  
  ## get common genes
  # removing genes with 0 counts everywhere
  gene1 = rownames(srat1@assays$RNA@counts)
  gene2 = rownames(srat2@assays$RNA@counts)
  common_genes = intersect(gene1,gene2)
  
  ## Calculate the fraction of total counts attributed to those genes that got removed
  srat1_unique_genes = rownames(srat1@assays$RNA@counts)[!rownames(srat1@assays$RNA@counts) %in% common_genes]
  srat1_frac_count_lost = sum(srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% srat1_unique_genes,])/sum(srat1@assays$RNA@counts)  
  
  srat2_unique_genes = rownames(srat2@assays$RNA@counts)[!rownames(srat2@assays$RNA@counts) %in% common_genes]
  srat2_frac_count_lost = sum(srat2@assays$RNA@counts[rownames(srat2@assays$RNA@counts) %in% srat2_unique_genes,])/sum(srat2@assays$RNA@counts)  
  
  message(sprintf('# genes in srat1: %d \n# genes in srat2: %d \n# genes in common: %d\n\nFrac count lost in srat1: %f \nFrac count lost in srat2: %f \n',
                  n_distinct(rownames(srat1@assays$RNA@counts)),
                  n_distinct(rownames(srat2@assays$RNA@counts)),
                  n_distinct(common_genes),
                  srat1_frac_count_lost,
                  srat2_frac_count_lost))
  
  
  
  cnt.mtx1 = srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% common_genes,]
  cnt.mtx2 = srat2@assays$RNA@counts[rownames(cnt.mtx1),]
  if(nrow(cnt.mtx1) != nrow(cnt.mtx2)){
    stop('There is issue with the merge!')
  }
  mtx = cbind(cnt.mtx1,cnt.mtx2)
  srat = CreateSeuratObject(mtx,meta.data = merged.mdat)
  
  return(srat)
}


### Converting annData to Seurat Obj
adata_to_srat = function(outdir,experiment_prefix){
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
  
  Convert(paste0(outdir, experiment_prefix, "_anndata.h5ad"),  
          dest = paste0(outdir, experiment_prefix, "_anndata.h5seurat"), overwrite = TRUE)
  seruat_object <- LoadH5Seurat(paste0(outdir, experiment_prefix, "_anndata.h5seurat"))
  
  return(seruat_object)
}



















