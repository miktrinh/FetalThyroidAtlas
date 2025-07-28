
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(edgeR)
library(singscore)
library(GenomicFeatures)
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)



# Combine TCGABiolinks::GDCprepare() and TCGABiolinks::readTranscriptomeProfiling() to make it such that it works with TCGA STAR-count transcriptomic data (instead of HTSeq-count)
process_TCGA_TranscriptomeProfiling = function(query,save = FALSE, save.filename, directory = "GDCdata", summarizedExperiment = TRUE,
                                               count_column_toKeep = c("unstranded","stranded_first","stranded_second","tpm_unstranded","fpkm_unstranded","fpkm_uq_unstranded")) {
  library(tidyverse)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  
  source <- ifelse(query$legacy, "legacy", "harmonized")
  files <- file.path(query$results[[1]]$project, source, gsub(" ", "_", query$results[[1]]$data_category), gsub(" ", "_", query$results[[1]]$data_type), gsub(" ", "_", query$results[[1]]$file_id), 
                     gsub(" ", "_", query$results[[1]]$file_name))
  files <- file.path(directory, files)
  
  cases <- ifelse(grepl("TCGA|TARGET", query$results[[1]]$project %>% 
                          unlist()), query$results[[1]]$cases, query$results[[1]]$sample.submitter_id)
  data.type = ifelse(!is.na(query$data.type), 
                     as.character(query$data.type), unique(query$results[[1]]$data_type))
  workflow.type = unique(query$results[[1]]$analysis_workflow_type)
  
  if (grepl("Gene Expression Quantification", data.type, ignore.case = TRUE)) {
    if (grepl("STAR", workflow.type)) {
      x <- plyr::alply(files, 1, function(f) {
        readr::read_tsv(file = f, col_names = TRUE, skip = 1,
                        progress = FALSE) %>%  
          dplyr::filter(!gene_id %in% c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous')) %>%
          dplyr::select(all_of(c('gene_id',count_column_toKeep)))
          #dplyr::rename_with(~ sampleName, -gene_id)
      }, .progress = "time")
      df <- x %>% purrr::reduce(left_join, by = "gene_id")
      if (!missing(cases)) 
        colnames(df)[-1] <- sapply(cases, function(x) {
          stringr::str_c(paste0(count_column_toKeep,'_'), x)
        }) %>% as.character()
      if (summarizedExperiment) 
        ## Adapted from TCGAbiolinks:::makeSEfromTranscriptomeProfilingSTAR()
        #df <- TCGAbiolinks:::makeSEfromTranscriptomeProfilingSTAR(df, cases, workflow.type)
        
        data = df
        size <- ncol(data)
        colData <- colDataPrepare(cases)
        gene.location <- get.GRCh.bioMart("hg38")
        gene.location <- gene.location[!duplicated(gene.location$ensembl_gene_id), 
                                       ]
        data$ensembl_gene_id <- as.character(gsub("\\.[0-9]*", "", 
                                                  data$gene_id))
        metrics <- subset(data, !grepl("ENSG", data$ensembl_gene_id))
        data <- subset(data, grepl("ENSG", data$ensembl_gene_id))
        found.genes <- table(data$ensembl_gene_id %in% gene.location$ensembl_gene_id)
        if ("FALSE" %in% names(found.genes)) 
          message(paste0("From the ", nrow(data), " genes we couldn't map ", 
                         found.genes[["FALSE"]]))
        data <- merge(data, gene.location, by = "ensembl_gene_id")
        assays <- list()
        for(count_type in count_column_toKeep){
          assays[[count_type]] <- data.matrix(data[, grep(paste0('^',count_type), colnames(data))])
        }
        
        assays <- lapply(assays, function(x) {
          colnames(x) <- NULL
          rownames(x) <- NULL
          return(x)
        })
        rowRanges <- GRanges(seqnames = paste0("chr", data$chromosome_name), 
                             ranges = IRanges(start = data$start_position, end = data$end_position), 
                             strand = data$strand, ensembl_gene_id = data$ensembl_gene_id, 
                             external_gene_name = data$external_gene_name, original_ensembl_gene_id = data$gene_id)
        names(rowRanges) <- as.character(data$ensembl_gene_id)
        rse <- SummarizedExperiment(assays = assays, rowRanges = rowRanges, 
                                    colData = colData)
        metadata(rse) <- metrics
        #return(rse)
        df = rse
    }
  }
  
  data = df
  
  
  if (summarizedExperiment & !is.data.frame(data)) {
    metadata(data) <- list(data_release = getGDCInfo()$data_release)
  }
  
  if ("samples" %in% colnames(data)) {
    if (any(duplicated(data$sample))) {
      message("Replicates found.")
      if (any(data$is_ffpe)) 
        message("FFPE should be removed. You can modify the data with the following command:\ndata <- data[,!data$is_ffpe]")
      print(as.data.frame(colData(data)[data$sample %in% 
                                          data$sample[duplicated(data$sample)], c("is_ffpe"), 
                                        drop = FALSE]))
    }
  }
  if (save) {
    if (missing(save.filename) & !missing(query)) 
      save.filename <- paste0(query$project, gsub(" ", 
                                                  "_", query$data.category), gsub(" ", "_", date()), 
                              ".RData")
    message(paste0("=> Saving file: ", save.filename))
    saveRDS(data, file = save.filename)
    message("=> File saved")
  }
  
  return(data)
}



is_named_vector <- function(x) {
  is.vector(x) && !is.null(names(x)) && all(names(x) != "")
}

compute_tpm <- function(count_matrix, gene_lengths) {
  # Ensure gene_lengths is in the same order as count_matrix rows
  gene_lengths <- gene_lengths[rownames(count_matrix)]
  
  # Convert lengths from bp to kilobases
  gene_lengths_kb <- gene_lengths / 1000
  
  # Step 1: Calculate RPK (Reads Per Kilobase)
  rpk <- sweep(count_matrix, 1, gene_lengths_kb, "/")
  
  # Step 2: Calculate per-sample scaling factor (sum of RPKs per sample)
  scaling_factors <- colSums(rpk)
  
  # Step 3: Calculate TPM
  tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6
  
  return(tpm)
}

#' @param bulk_sources: named vector of file paths to the RDS objects of each dataset
import_bulkRNA_thyroid = function(bulk_sources = c('He2021'='Data/published_bulkRNAseq/He_etal_21/aPTC_He_2021_se.RDS','Lee2021'=NULL,'Yoo2016'=NULL,
                                                   'Lee2024' = 'Data/published_bulkRNAseq/Lee_etal_24/aPTC_Lee_2024_se.RDS',
                                                   'TCGA_Thyroid'='Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS','snPaedThyroid'=NULL,'scaThy'=NULL,'GTEx_Thyroid'=NULL,'StJudes_Thyroid'=NULL,'fAdr'=NULL,
                                                   'inhouse'='Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS'),
                                  gene_map=gene_map){
  library(readxl)
  library(SummarizedExperiment)
  library(edgeR)
  library(singscore)
  library(GenomicFeatures)
  
  
  
  if (!is.vector(bulk_sources) || is.null(names(bulk_sources)) || any(names(bulk_sources) == "")) {
    stop("bulk_sources must be a named vector.")
  }
  # Remove datasets with no file path provided
  bulk_sources = bulk_sources[!is.null(bulk_sources)]
  
  
  ##----- 1. Import bulk counts -----##
  bulk_samples = dplyr::tibble()
  sce_list = list()
  raw_count = data.frame()
  tpm_count = data.frame()
  
  
  
  if('He2021' %in% names(bulk_sources)){
    se_path = bulk_sources[['He2021']]
     
    if(!file.exists(se_path)){
      rawCnt = read.delim('Data/published_bulkRNAseq/He_etal_21/GSE165724_counts_74Samples.tsv.gz',sep = '\t')
      
      # Get sample metdata
      gse <- GEOquery::getGEO("GSE165724", GSEMatrix = TRUE)
      # If multiple platforms, choose the first one
      gse <- gse[[1]]
      sample_metadata <- pData(gse)
      
      colnames(rawCnt) <- c('Geneid',sample_metadata$geo_accession)
      
      sample_metadata$sampleID = sample_metadata$geo_accession
      sample_metadata$source = 'He_2021'
      sample_metadata$sampleName = sample_metadata$geo_accession
      sample_metadata$cancerType = ifelse(grepl('tissue type: Normal$',sample_metadata$characteristics_ch1.2),'Normal',
                                          ifelse(grepl('tissue type: Normal ',sample_metadata$characteristics_ch1.2),'Normal.adj',
                                            ifelse(grepl('tissue type: PTC',sample_metadata$characteristics_ch1.2),'PTC','Others')))
      sample_metadata$age = as.numeric(gsub('age: ','',sample_metadata$characteristics_ch1.4))
      sample_metadata$sex = gsub('gender: ','',sample_metadata$characteristics_ch1.3)
      
      row_data = data.frame('Gene_id' = rawCnt[,'Geneid'])    
      d = as.data.frame(gene_map[gene_map$gene_id %in% gsub('\\.\\d+$','',rawCnt$Geneid)])
      d$gene_id_stripped = d$gene_id
      row_data = row_data %>%
        mutate(gene_id_stripped = gsub('\\.\\d+$','',row_data$Gene_id)) %>% 
        left_join(d, by = "gene_id_stripped")
      rownames(row_data) = row_data$gene_id_stripped
      
      rownames(rawCnt) = row_data$gene_id_stripped
      
      # Import gene length
      geneLength = readRDS('Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS') %>% rowData()
      geneLength = geneLength[(geneLength$ensemblID %in% rownames(rawCnt)) & 
                                geneLength$effLength >0,]
      
      rawCnt = rawCnt[rownames(rawCnt) %in% geneLength$ensemblID,colnames(rawCnt) != 'Geneid']
      row_data = row_data[row_data$gene_id_stripped %in% geneLength$ensemblID,]
      row_data = cbind(row_data,geneLength[match(row_data$gene_id_stripped,geneLength$ensemblID),])
        
        
      geneLength = geneLength[match(rownames(rawCnt),geneLength$ensemblID),]
      geneLength = setNames(geneLength$effLength,geneLength$ensemblID)
      
      tpmCnt = compute_tpm(count_matrix=rawCnt, gene_lengths=geneLength) 
      
      rse <- SummarizedExperiment(assays = list(counts_raw = rawCnt,
                                                counts_tpm = tpmCnt), 
                                  rowData = row_data[match(rownames(rawCnt),row_data$gene_id_stripped),], 
                                  colData = sample_metadata)
      saveRDS(rse,se_path)
    }
    
    he21_se = readRDS(se_path)
    sce_list[['He2021']] = he21_se
    he21_mdat = as.data.frame(colData(he21_se))
    
    bulk_samples = rbind(bulk_samples,he21_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    he21_rawCnt = assays(he21_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = he21_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(he21_rawCnt))
      raw_count = cbind(raw_count[genesToKeep,],he21_rawCnt[genesToKeep,])
    }
  }
  
  
  if('Lee2024' %in% names(bulk_sources)){
    se_path = bulk_sources[['Lee2024']]
    
    if(!file.exists(se_path)){
      files = list.files('Data/published_bulkRNAseq/Lee_etal_24/',pattern = 'ReadsPerGene.out.tab.gz$',full.names = T,recursive = F)
      
      x <- plyr::alply(files, 1, function(f) {
        sampleName = gsub('_.*$','',basename(f))
        readr::read_tsv(file = f, col_names = F, skip = 4,progress = FALSE) %>%  
          magrittr::set_names(c("gene_id", "unstranded", "stranded_forward", "stranded_reverse")) %>%
          dplyr::select(all_of(c('gene_id',"unstranded"))) %>% 
          dplyr::rename_with(~ sampleName, -gene_id)
      }, .progress = "time")
      rawCnt <- x %>% purrr::reduce(left_join, by = "gene_id")
      rawCnt <- column_to_rownames(rawCnt,'gene_id')
      
      # Get sample metdata
      Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
      gse <- GEOquery::getGEO("GSE213647", GSEMatrix = TRUE)
      length(gse)
      # If multiple platforms, choose the first one
      gse <- gse[[1]]
      sample_metadata <- pData(gse)
      
      sample_metadata$sampleID = sample_metadata$geo_accession
      sample_metadata$source = 'Lee_2024'
      sample_metadata$sampleName = sample_metadata$geo_accession
      sample_metadata$cancerType = ifelse(sample_metadata$`cell type:ch1` == 'Normal','Normal',
                                          ifelse(sample_metadata$`cell type:ch1` == 'Tumor',sample_metadata$`cell subtype:ch1`,'Others'))
      sample_metadata$age = 'adult'
      sample_metadata$sex = '-'
      
      
      row_data = data.frame('gene_id' = rownames(rawCnt))    
      d = as.data.frame(gene_map[gene_map$gene_id %in% gsub('\\.\\d+$','',rownames(rawCnt))])
      d$gene_id_stripped = d$gene_id
      row_data = row_data %>%
        mutate(gene_id_stripped = gsub('\\.\\d+$','',row_data$gene_id)) %>% 
        left_join(d[,colnames(d)!='gene_id'], by = "gene_id_stripped")
      rownames(row_data) = row_data$gene_id_stripped
      
      rownames(rawCnt) = row_data$gene_id_stripped
      
      # Import gene length
      geneLength = readRDS('Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS') %>% rowData()
      geneLength = geneLength[(geneLength$ensemblID %in% rownames(rawCnt)) & 
                                geneLength$effLength >0,]
      
      rawCnt = rawCnt[rownames(rawCnt) %in% geneLength$ensemblID,]
      row_data = row_data[row_data$gene_id_stripped %in% geneLength$ensemblID,]
      row_data = cbind(row_data,geneLength[match(row_data$gene_id_stripped,geneLength$ensemblID),])
      
      
      geneLength = geneLength[match(rownames(rawCnt),geneLength$ensemblID),]
      geneLength = setNames(geneLength$effLength,geneLength$ensemblID)
      
      tpmCnt = compute_tpm(count_matrix=rawCnt, gene_lengths=geneLength) 
      
      rse <- SummarizedExperiment(assays = list(counts_raw = rawCnt,
                                                counts_tpm = tpmCnt), 
                                  rowData = row_data[match(rownames(rawCnt),row_data$gene_id_stripped),], 
                                  colData = sample_metadata)
      
      saveRDS(rse,se_path)
    }
    
    lee24_se = readRDS(se_path)
    samples_toKeep = colData(lee24_se)
    samples_toKeep = samples_toKeep$sampleID[samples_toKeep$`cell subtype:ch1`=='PTC']
    lee24_se = lee24_se[,samples_toKeep]
    sce_list[['Lee2024']] = lee24_se
    lee24_mdat = as.data.frame(colData(lee24_se))
    
    bulk_samples = rbind(bulk_samples,lee24_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    lee24_rawCnt = assays(lee24_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = lee24_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(lee24_rawCnt))
      raw_count = cbind(raw_count[genesToKeep,],lee24_rawCnt[genesToKeep,])
    }
  }
    
    
    
  if('Lee2021' %in% names(bulk_sources)){
    se_path = '~/lustre_mt22/Thyroid/Data/Lee_etal_2021/pedPTC_Lee_2021_se.RDS'
    pedPTC_se = readRDS(se_path)
    sce_list[['Lee2021']] = pedPTC_se
    pedPTC_mdat = as.data.frame(colData(pedPTC_se))
    
    pedPTC_mdat$sampleName = pedPTC_mdat$Sample.Name
    pedPTC_mdat$age = pedPTC_mdat$Age
    bulk_samples = rbind(bulk_samples,pedPTC_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    pedPTC_rawCnt = assays(pedPTC_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = pedPTC_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(pedPTC_rawCnt))
      raw_count = cbind(raw_count[genesToKeep,],pedPTC_rawCnt[genesToKeep,])
    }
    
  }
  
  
  if('Yoo2016' %in% names(bulk_sources)){
    se_path = bulk_sources[['Yoo2016']]
    
    
    
    aThy_se = readRDS(se_path)
    sce_list[['Yoo2016']] = aThy_se
    
    aThy_mdat = as.data.frame(colData(aThy_se))
    bulk_samples = rbind(bulk_samples,aThy_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    aThy_rawCnt = assays(aThy_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = aThy_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(aThy_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],aThy_rawCnt[genesToKeep,])
    }
    
  }
  
  if('TCGA_Thyroid' %in% names(bulk_sources)){
    #tcga_se_path_m2 = '~/lustre_mt22/Thyroid/Data/TCGA_Thyroid/TCGA_Thyroid_gdc0923_sce.RDS'
    tcga_se_path = bulk_sources['TCGA_Thyroid']
    
    
    
    # TCGAbiolinks docs: https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
    if(!file.exists(tcga_se_path)){
      library(TCGAbiolinks)
      
      query <- TCGAbiolinks::GDCquery(
        project = "TCGA-THCA",  
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
      )
      
      GDCdownload(query,directory = 'Data/published_bulkRNAseq/TCGA_Thyroid/GDCdata')
      #data <- GDCprepare(query,directory = "Data/published_bulkRNAseq/TCGA_Thyroid/GDCdata", summarizedExperiment = FALSE)
      data <- process_TCGA_TranscriptomeProfiling(query,save = TRUE, save.filename=tcga_se_path, directory = "Data/published_bulkRNAseq/TCGA_Thyroid/GDCdata", summarizedExperiment = TRUE,
                                                     count_column_toKeep = c("unstranded","stranded_first","stranded_second","tpm_unstranded","fpkm_unstranded","fpkm_uq_unstranded"))
    }
    
    tcga_se = readRDS(tcga_se_path)
    tcga_mdat = as.data.frame(colData(tcga_se))
    # Keep only samples of the relevant cancer types only
    # According to this paper here: https://pmc.ncbi.nlm.nih.gov/articles/PMC10298340/
    samples_toKeep = rownames(tcga_mdat[grepl('Papillary carcinoma|Papillary adenocarcinoma',tcga_mdat$primary_diagnosis) &
                                 tcga_mdat$tissue_or_organ_of_origin %in% c('Thyroid gland',"Lymph node, NOS"),])
    
    tcga_se = tcga_se[,samples_toKeep]
    sce_list[['TCGA_Thyroid']] = tcga_se
    assays(sce_list[['TCGA_Thyroid']])[['counts_tpm']] = assays(sce_list[['TCGA_Thyroid']])[['tpm_unstranded']]
    
    tcga_mdat = as.data.frame(colData(tcga_se))
    tcga_mdat$sampleID= rownames(tcga_mdat)
    tcga_mdat$source = 'TCGA_Thyroid'
    tcga_mdat$sampleName = rownames(tcga_mdat)
    tcga_mdat$cancerType = ifelse(tcga_mdat$tissue_type == 'Normal','Normal',
                                  paste0('PTC_',tcga_mdat$classification_of_tumor))
    tcga_mdat$age = tcga_mdat$age_at_diagnosis
    tcga_mdat$sex = tcga_mdat$gender
    
    bulk_samples = rbind(bulk_samples,tcga_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    
    
    tcga_rawCnt = assays(tcga_se)[['unstranded']]
    if(nrow(raw_count) == 0){
      raw_count = tcga_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(tcga_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],tcga_rawCnt[genesToKeep,])
    }
    
  }
  
  
  if('snPaedThyroid' %in% names(bulk_sources)){
    se_path = '~/lustre_mt22/Thyroid/Data/paedThyroid_snRNAseq_pseudobulk_sce.RDS'
    pThyroid_sce = readRDS(se_path)
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
  
  if('scfThy' %in% names(bulk_sources)){
    se_path = '~/lustre_mt22/Thyroid/Data/fThy_scRNAseq_pseudobulk_sce.RDS'
    scfThy_sce = readRDS(se_path)
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
  
  
  
  if('scaThy' %in% names(bulk_sources)){
    se_path = '~/lustre_mt22/Thyroid/Data/aThy_scRNAseq_pseudobulk_sce.RDS'
    scaThy_sce = readRDS(se_path)
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
  
  
  if('GTEx_Thyroid' %in% names(bulk_sources)){
    gtex_se_path = '~/lustre_mt22/Thyroid/Data/GTEx_Thyroid/GTEx_Thyroid_se.RDS'
    gtex_se = readRDS(gtex_se_path)
    sce_list[['GTEx_Thyroid']] = gtex_se
    
    gtex_mdat = as.data.frame(colData(gtex_se))
    bulk_samples = rbind(bulk_samples,gtex_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    gtex_rawCnt = assays(gtex_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = gtex_rawCnt
    }else{
      rowDat = rowData(gtex_se)
      #rownames(gtex_rawCnt) = rowDat$ensID[match(rownames(gtex_rawCnt),rowDat$gene_id)]
      genesToKeep = intersect(rownames(raw_count),rownames(gtex_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],gtex_rawCnt[genesToKeep,])
    }
  }
  
  
  if('StJudes_Thyroid' %in% names(bulk_sources)){
    stJudes_se_path = '~/lustre_mt22/Thyroid/Data/StJudes_Thyroid/StJudes_Thyroid_230921_se.RDS'
    stJudes_se = readRDS(stJudes_se_path)
    sce_list[['StJudes_Thyroid']] = stJudes_se
    
    stJudes_mdat = as.data.frame(colData(stJudes_se))
    bulk_samples = rbind(bulk_samples,stJudes_mdat[,c('sampleID','source','sampleName','cancerType','age','sex')])
    
    stJudes_rawCnt = assays(stJudes_se)[['counts_raw']]
    if(nrow(raw_count) == 0){
      raw_count = stJudes_rawCnt
    }else{
      genesToKeep = intersect(rownames(raw_count),rownames(stJudes_rawCnt))
      print(length(genesToKeep))
      raw_count = cbind(raw_count[genesToKeep,],stJudes_rawCnt[genesToKeep,])
    }
  }
  
  if('inhouse' %in% names(bulk_sources)){
    se_path = bulk_sources['inhouse']
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
  
  if('fAdr' %in% names(bulk_sources)){
    se_path = '~/lustre_mt22/Thyroid/Data/fAdr.SCPs_scRNAseq_pseudobulk_sce.RDS'
    scfAdr_sce = readRDS(se_path)
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
  tpm_count = data.frame()
  for(i in 1:length(sce_list)){
    print(i)
    tpm = assays(sce_list[[i]])[['counts_tpm']]
    # if(names(sce_list)[i] == 'GTEx_Thyroid'){
    #   rowDat = rowData(sce_list[[i]])
    #   rownames(tpm) = rowDat$ensID[match(rownames(tpm),rowDat$gene_id)]
    # }
    
    if(nrow(tpm_count) == 0){
      tpm_count = tpm
    }else{
      genesToKeep = intersect(rownames(tpm_count),rownames(tpm))
      print(length(genesToKeep))
      tpm_count = cbind(tpm_count[genesToKeep,],tpm[genesToKeep,])
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
  bulk_dge = DGEList(counts = raw_count, genes = rownames(raw_count))
  
  ## Plot library size
  libSize = colSums(raw_count)
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
  
  raw_count_geneDetection = (raw_count>0)
  raw_count_geneDetection = colSums(raw_count_geneDetection)
  bulk_samples$nGene = raw_count_geneDetection[match(bulk_samples$sampleID,names(raw_count_geneDetection))]
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
  
  raw_count = raw_count[keep, ]
  bulk_dge = DGEList(counts = raw_count, genes = rownames(raw_count))
  tpm_count = tpm_count[rownames(tpm_count) %in% rownames(raw_count),]
  
  cpmCnt = edgeR::cpm(bulk_dge, log = TRUE, prior.count = 1)
  # lapply(split(colnames(cpmCnt),bulk_samples$source),function(e){
  #   dataset = unique(bulk_samples$source[bulk_samples$sampleID %in% e])
  #   print(hist(cpmCnt[,e],main = dataset))})
  
  return(list('raw_count'=raw_count,
              'tpm_count'=tpm_count,
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

