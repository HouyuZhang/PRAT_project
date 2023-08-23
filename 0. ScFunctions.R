#===============================================================================
# This script stores functions to analyze sc/snRNA-seq dataset

# Version 1.4
# Created by Houyu Zhang
# Issue report on houyuzhang@stu.pku.edu.cn
# Copyright (c) 2023 __CarlosLab@PKU__. All rights reserved.
#===============================================================================
#Global setup
R.utils::setOption("clusterProfiler.download.method",'auto')
#General packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(igraph))
library("data.table")
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsankey))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(EnhancedVolcano))
#Single cell packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(SeuratDisk)) #remotes::install_github("mojaveazure/seurat-disk")
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(monocle3)) #devtools::install_github('cole-trapnell-lab/monocle3')
suppressPackageStartupMessages(library(DoubletFinder)) #remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(CellChat)) #devtools::install_github("sqjin/CellChat") also, circlize
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(clusterProfiler)) #BiocManager::install("clusterProfiler")
suppressPackageStartupMessages(library(ReactomePA)) #BiocManager::install("ReactomePA")
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(org.Mm.eg.db)) #BiocManager::install("org.Mm.eg.db")
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(wpa))

#===============================================================================
# Pre-set color schemes for snRNA-seq datasets
#===============================================================================
hughie_color <- c("#0072B5FF","#BC3C29FF", "#E18727FF", "#20854EFF", "#6F99ADFF", "#7876B1FF","#EE4C97FF", #nemj scheme
                  "#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", #npg scheme
                  "#3B4992FF", "#631879FF", "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF") #AAAS scheme
PairedColor <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")


#===============================================================================
# Custome functions for sc/snRNA-seq analysis 
#===============================================================================
##' A warpper to run Seurat preprocessing
##' @param Seurat_Obj A Seurat object
##' @param NormalizationMethod
##' @param removeDoublet Logic value to act on doublet remove
##' @param DoubletPercent Percentage of doublet in your dataset (8% for 20K cells loaded in 10X)
##' @param species Specify species(mouse, rat, human) you researched (default: mouse)
##' @param Filter_nFeature_RNA threshold of minimal and maximal genes detected in each cell
##' @param Filter_nCount_RNA threshold of minimal and maximal UMI counts detected in each cell
##' @param MitoPercent threshold to remove cells based on mitochondria reads (default: 15%)
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' RunSeurat(ObjList = c(rATM1mo_Obj, rATM2mo_Obj, rATM6mo_Obj), OutPrefix = rATM_prefix, NormalizationMethod = "VST",
##'           Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)

PlotFeature <- function(ObjList, Prefix){
  pdf(paste0(Prefix,".pdf"), width = 9, height = 6)
  for(i in 1:length(ObjList)){p <- VlnPlot(ObjList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3); plot(p)}
  dev.off()
}

RunSeurat <- function(ObjList, OutPrefix, NormalizationMethod = "VST", removeDoublet = T, species = "mouse",
                      Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15, DoubletPercent = 0.08){

  if(species == "mouse"){mtPatterm <- "^mt-"
  } else if(species == "rat"){mtPatterm <- "^Mt-"
  } else if(species == "human"){mtPatterm <- "^MT-"
  }
  
  # Plot data metrics before filtering
  PlotFeature(ObjList, Prefix = paste0(OutPrefix,"BeforeFilter_",Sys.Date()))
  # future::plan("multisession", workers = 20) # do parallel

  if(NormalizationMethod == "VST"){
    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      #1. Low-quality cells or empty droplets will often have very few genes. 
      #2. Cell doublets or multiplets may exhibit an aberrantly high gene count. 
      #3. Low-quality / dying cells often exhibit extensive mitochondrial contamination
      Obj[["percent.mt"]] <- PercentageFeatureSet(Obj, pattern = mtPatterm)
      Obj <- subset(Obj, subset = nFeature_RNA > Filter_nFeature_RNA[1] & nFeature_RNA < Filter_nFeature_RNA[2] 
                    & nCount_RNA > Filter_nCount_RNA[1] & nCount_RNA < Filter_nCount_RNA[2] & percent.mt < MitoPercent)
      Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
      Obj <- ScaleData(Obj) %>% RunPCA(dims = 1:15) %>% RunUMAP(dims = 1:15) %>% 
        FindNeighbors(reduction = "pca", dims = 1:15) %>% FindClusters(resolution = seq(from=0, by=0.05, length=3))
    })
    
    pdf(paste0(OutPrefix,"UMAP-ByDoublets_",Sys.Date(),".pdf"), width = 8, height = 8)
    if(removeDoublet){
      for (i in 1:length(ObjList)){
        ## pK Identification
        sweep.res.list <- paramSweep_v3(ObjList[[i]], PCs = 1:15, sct = F)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal.pk <- as.numeric(as.vector(bcmvn[which(bcmvn$BCmetric == max(bcmvn$BCmetric)),]$pK))
        
        ## Homotypic Doublet Proportion Estimate
        annotations <- ObjList[[i]]@meta.data$RNA_snn_res.0.1
        homotypic.prop <- modelHomotypic(annotations)
        nExp_poi <- round(DoubletPercent*nrow(ObjList[[i]]@meta.data)) ## Assuming 8% doublet formation rate based 10X guideline when load 20,000 cells
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## Run DoubletFinder with varying classification stringencies
        ObjList[[i]] <- doubletFinder_v3(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
        # ObjList[[i]] <- doubletFinder_v3(ObjList[[i]], PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",optimal.pk,"_",nExp_poi), sct = FALSE)
        Idents(ObjList[[i]]) <- ObjList[[i]]@meta.data[,paste0("DF.classifications_0.25_",optimal.pk,"_",nExp_poi)]
        p <- DimPlot(ObjList[[i]], reduction = "umap", label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white");plot(p)
        Idents(ObjList[[i]]) <- ObjList[[i]]$seurat_clusters
        toRemove <- rownames(ObjList[[i]]@meta.data[ObjList[[i]]@meta.data[,ncol(ObjList[[i]]@meta.data)] == "Doublet",])
        ObjList[[i]] <- ObjList[[i]][,!colnames(ObjList[[i]]) %in% toRemove]
        cat(paste0(length(toRemove)," doublets were removed !!!\n"))
      }
    }
    dev.off()

    ObjList <- lapply(X = ObjList, FUN = function(Obj) {
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
      Obj <- ScaleData(Obj)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
    anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
    combined <- IntegrateData(anchorset = anchors)
    DefaultAssay(combined) <- "integrated"
  }
  
  # Plot data metrics after filtering
  PlotFeature(ObjList, Prefix = paste0(OutPrefix,"AfterFilter_",Sys.Date()))

  if(NormalizationMethod == "SCT"){pc = 15}else{pc = 30}
  combined <- ScaleData(combined) %>% RunPCA(npcs = pc) %>% RunUMAP(reduction = "pca", dims = 1:pc) %>% 
    FindNeighbors(reduction = "pca", dims = 1:pc) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(combined)
  ggsave(paste0(OutPrefix,"clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(combined) <- combined$integrated_snn_res.0.1
  pdf(paste0(OutPrefix,"AllFeatureQC_",Sys.Date(),".pdf"), width = 12, height = 12)
  p <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1); plot(p); dev.off()
  saveRDS(combined, file = paste0(OutPrefix,"Seurat.rds"))
  
  # Visualization of the clusters landscape
  DimPlot(combined, reduction = "umap", group.by = c("ident","orig.ident"), label = T, repel = T, label.box =T, cols = hughie_color, label.color = "white")
  ggsave(paste0(OutPrefix,"UMAP-ByClusters_",Sys.Date(),".pdf"), width = 16, height = 8)
  
  DimPlot(combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = T, repel = T, cols = hughie_color)
  ggsave(paste0(OutPrefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 12)
  
  Plot_Cell_compoistion(Seurat_Obj = combined, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = OutPrefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = combined, ColorUse = hughie_color, OutPrefix = OutPrefix)
  
  # future::plan("multisession", workers = 4)
  return(combined)
}

##' Function to caculate cell type module score based on marker genes
##' @param Seurat_Obj A Seurat object
##' @param Seruat_DEG_file DEG file for reference cell types, returned by the FindAllMarkers function
##' @param NtopGene Number of top genes (ranked by log2 FC) used as a module (Default: 50)
##' @param Ylimits Range for Y-axis on plotting
##' @param convertGeneName convert gene names when comparing different species, H2M indicate human to mouse, etc
##' @param ColorUse Color vector 
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Calculate_moduleScore(Seurat_Obj = subset(iBAT_Ad_combined, subset = orig.ident == "iBAT-2mo"), ColorUse = iBAT_Ad_colors, NtopGene = 50, Ylimits = c(0, 1),
##'                       Seruat_DEG_file = "../1. rATM_126mo/rATM_Adipocytes_2mo_Allmarkers_2023-07-05.csv", OutPrefix = paste0(iBAT_Ad_prefix, "2mo_mPRAT_"))
Calculate_moduleScore <- function(Seurat_Obj, Seruat_DEG_file, NtopGene = 50, OutPrefix, Ylimits = c(-2,2), convertGeneName = FALSE, ColorUse = hughie_color){
  DEG_matrix <- read.table(Seruat_DEG_file, sep = ",", header = T, row.names = 1)
  DEG_matrix1 <- DEG_matrix %>% group_by(cluster) %>% slice_max(n = NtopGene, order_by = avg_log2FC)
  DEG_list <- split(DEG_matrix1$gene, f = DEG_matrix1$cluster)
  
  if(convertGeneName == "H2M"){DEG_list <- lapply(DEG_list, HumanGeneToMouseGene)}
  if(convertGeneName == "M2H"){DEG_list <- lapply(DEG_list, MouseGeneToHumanGene)}
  
  pdf(paste0(OutPrefix,"UMAP-ModuleScore_",NtopGene,"TopMarkers_",Sys.Date(),".pdf"), width = 8, height = 10)
  for(CellTypes in names(DEG_list)){
    Seurat_Obj <- AddModuleScore_UCell(Seurat_Obj, features = list(DEG_list[[CellTypes]]), name = paste0(make.names(CellTypes),"_score"), assay = "RNA")
    p <- FeaturePlot(Seurat_Obj, features = paste0("signature_1",make.names(CellTypes),"_score"), min.cutoff = "q10", order = T, pt.size=0.3) + 
      scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
    plot(p)
  }
  p1 <- VlnPlot(Seurat_Obj, features = paste0("signature_1",make.names(names(DEG_list)),"_score"), split.by = "orig.ident", stack = T, flip = T,
                fill.by = "ident", cols = hughie_color) + FontSize(x.text = 16, y.text = 16) + scale_y_continuous(limits = Ylimits)
  p2 <- VlnPlot(Seurat_Obj, features = paste0("signature_1",make.names(names(DEG_list)),"_score"), stack = T, flip = T,
                fill.by = "ident", cols = hughie_color) + FontSize(x.text = 16, y.text = 16) + scale_y_continuous(limits = Ylimits)
  plot(p1); plot(p2)
  dev.off()
  
  Seurat_Obj$activeId <- Idents(Seurat_Obj)
  sigMat <- Seurat_Obj@meta.data[,c("activeId",paste0("signature_1",make.names(names(DEG_list)),"_score"))]
  sigMat <- reshape2::melt(sigMat,c("activeId"))
  sigMat$variable <- gsub("signature_1","",sigMat$variable)
  
  p3 <- ggplot(sigMat, aes(x = activeId, y = value)) + 
    geom_violin(aes(fill = activeId),width=1, linewidth =0.1) +
    geom_boxplot(color="black", alpha=0.2, width=0.1, linewidth =0.1, outlier.shape = NA) +
    stat_compare_means(method = "t.test", comparisons = combn(as.vector(unique(sigMat$activeId)), 2, simplify = FALSE)) +
    # stat_compare_means(label.y = 0.6) +
    facet_wrap(~variable, ncol=6, scales = "free_x") +
    scale_fill_manual(values = ColorUse) +
    expand_limits(y=0) +
    labs(x = "", y = "Module score") +
    theme_bw() +
    theme(
      plot.title = element_text(color="black", size=10, face="bold"),
      axis.title.x = element_text(color="black", size=10, face="bold"),
      axis.text.x = element_text(color="black", size=10, face="bold", angle = 45, hjust = 1),
      axis.title.y = element_text(color="black", size=10, face="bold"),
      axis.text.y = element_text(color="black", size=10, face="bold"),
      legend.title = element_text(color="black", size=10, face="bold"), 
      legend.text = element_text(color="black", size=10, face="bold"),
      strip.text.x = element_text(size = 8, face="bold")
    )
  plot(p3)
  ggsave(paste0(OutPrefix,"UMAP-ModuleScore_",NtopGene,"TopMarkers_Integrated",Sys.Date(),".pdf"), width = 14, height = 10)

  return(Seurat_Obj)
}

##' A wrapper to do DEG and plotting
##' @param Seurat_Obj A Seurat object
##' @param NMarkers Number of markers for violin and heatmap
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' MarkersPlot(Seurat_Obj = iWAT_combined, ColorUse = iWAT_HFD_color_maps, OutPrefix = "_iWAT_combined_",NMarkers = 20)
MarkersPlot <- function(Seurat_Obj, OutPrefix, NMarkers = 20, ColorUse = hughie_color){
  
  markers <- FindAllMarkers(Seurat_Obj, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25, assay = "RNA") %>% filter(p_val_adj < 0.05) %>% filter(!str_detect(gene,"^mt-"))
  write.csv(markers, file = paste0(OutPrefix, "Allmarkers_",Sys.Date(),".csv"), quote = F)
  TopMarkersList <- markers %>% group_by(cluster) %>% slice_max(n = NMarkers, order_by = avg_log2FC) %>% 
    arrange(match(cluster, levels(Idents(Seurat_Obj)))) %>% as.data.frame()
  
  # Plot all markers
  VlnPlot(Seurat_Obj, features = unique(TopMarkersList$gene), stack=T, flip=T, cols = ColorUse, fill.by = "ident") + NoLegend()
  ggsave(paste0(OutPrefix,"Allmarkers_Vlnplot_",Sys.Date(),".pdf"), width = 12, height = 45)
  
  # DotPlot(Seurat_Obj, features = unique(TopMarkersList$gene), split.by = "orig.ident", cols = hughie_color, scale = T) + RotatedAxis()
  # ggsave(paste0(OutPrefix,  "Allmarkers_Dotplot_",Sys.Date(),".pdf"), width = 16, height = 10)
  
  DoHeatmap(Seurat_Obj, features = unique(TopMarkersList$gene), group.colors = ColorUse)
  ggsave(paste0(OutPrefix, "Allmarkers_Heatmap_",Sys.Date(),".pdf"), width = 16, height = 20)
  return(TopMarkersList)
}

##' Function to plot cell quality metrice from a Seurat Obj
##' @param Seurat_Obj A Seurat object
##' @param ColorUse Color vector 
##' @param OutPrefix Prefix for output file
PlotObjMetrices <- function(Seurat_Obj, OutPrefix, ColorUse = hughie_color){
  met <- Seurat_Obj@meta.data[c("orig.ident","nCount_RNA","percent.mt","nFeature_RNA")] %>% reshape2::melt(c("orig.ident"))
  met$variable <- factor(met$variable, levels = c("nCount_RNA","nFeature_RNA","percent.mt" ))
  
  p <- ggplot(met, aes(x = orig.ident, y = value)) + 
    geom_violin(aes(fill = orig.ident), width=1, linewidth=0.1) +
    geom_boxplot(color="black", alpha=0.1, width=0.1, outlier.colour = NA, linewidth=0.1) +
    facet_wrap(~variable, ncol=6, scales = "free") +
    # facetted_pos_scales(y = list(variable == "nCount_RNA" ~ scale_y_continuous(limits = c(0, 25000),  breaks = seq(0, 25000, 5000)),
    #                              variable == "nFeature_RNA" ~ scale_y_continuous(limits = c(0, 5000),  breaks = seq(0, 5000, 1000)),
    #                              variable == "percent.mt" ~ scale_y_continuous(limits = c(0, 10),  breaks = seq(0, 10, 2)))) +
    scale_fill_manual(values = ColorUse) +
    labs(x = "", y = "") + theme_bw() +
    theme(
      axis.title.x = element_text(color="black", size=8, face="bold"),
      axis.text.x = element_text(color="black", size=8, face="bold", angle = 45, hjust = 1),
      axis.title.y = element_text(color="black", size=8, face="bold"),
      axis.text.y = element_text(color="black", size=8, face="bold"),
      strip.text.x = element_text(size = 8, face="bold")
    )
  plot(p)
  ggsave(paste0(OutPrefix,"snRNA-seq_metrics.pdf"), width = 8, height = 4)
}

##' Read h5 files
##' @param filename h5 filename
##' @examples
##'
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    message('CellRanger version 3+ format H5')
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    message('CellRanger version 2 format H5')
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      repr="T"
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'CsparseMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

##' Print cell composition
GetCellComposition <- function(Seurat_Obj){
  CompositionFreq <- table(Seurat_Obj$orig.ident, Seurat_Obj@active.ident) %>% as.data.frame() 
  mat <- CompositionFreq %>% group_by(Var1) %>% mutate(Freq = Freq/sum(Freq) * 100) 
  mat$Var2 <- factor(mat$Var2, levels = rev(levels(mat$Var2)))
  colnames(mat) <- c("Sample","CellType","Perc(%)")
  print(as.data.frame(mat),digits=2)
  # gridExtra::grid.table(mat)
}

##' Plot multiple marker genes expression
##' @param Seurat_Obj A Seurat object
##' @param ColorUse1 List of colors for cell types
##' @param ColorUse2 List of colors for sample types
##' @param OutPrefix Prefix for output file
##' @examples
##'
##' Plot_Cell_compoistion(Seurat_Obj = iWAT_combined, OutPrefix = "_iWAT_combined_", ColorUse = iWAT_HFD_color_maps)
Plot_Cell_compoistion <- function(Seurat_Obj, ColorUse1, ColorUse2, OutPrefix){
  
  CompositionFreq <- table(Seurat_Obj$orig.ident, Seurat_Obj@active.ident) %>% as.data.frame() 
  
  mat <- CompositionFreq %>% group_by(Var1) %>% mutate(Freq = Freq/sum(Freq) * 100) 
  mat$Var2 <- factor(mat$Var2, levels = rev(levels(mat$Var2)))
  p2 <- ggplot(mat,aes(x=Var1, y=Freq, fill = Var2, label = paste0(round(Freq, 1), "%"))) + 
    geom_bar(stat="identity", width = 0.7, color = "black") +
    geom_text(position=position_stack(vjust=0.5), color="white") + 
    # scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = ColorUse1) + 
    coord_flip() +  scale_x_discrete(limits = rev(levels(mat$Var1))) + theme_bw() + 
    labs(x="", y="Percentage of each cell type") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x = element_text(color="black", size=16, face="bold"),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=16, face="bold")
    )
  plot(p2)
  ggsave(paste0(OutPrefix, "cellNumbersPerCluster_filled_",Sys.Date(),".pdf"), width = 10, height = 5)
  
  p4<-ggplot(mat,aes(x="", y=Freq, fill = Var2, label = paste0(round(Freq, 1), "%"))) + 
    geom_bar(stat="identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    geom_text(position=position_stack(vjust=0.5), color="white",size=8) + 
    scale_fill_manual(values = ColorUse1) + 
    facet_wrap(~Var1) +
    theme_void() + 
    labs(x="", y="Percentage of each cell type") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x=element_blank(),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=16, face="bold")
    )
  plot(p4)
  ggsave(paste0(OutPrefix, "cellNumbersPerCluster_Pie_",Sys.Date(),".pdf"), width = 8, height = 4)
  
  NCellType <- length(unique(Seurat_Obj@active.ident))
  Quan <- table(Seurat_Obj$orig.ident,Seurat_Obj@active.ident) %>% as.data.frame.matrix()
  Quan$Sum <- rowSums(Quan)
  Quan <- Quan*100/Quan$Sum
  Quan <- Quan[,-ncol(Quan)]
  Quan$Samples <- rownames(Quan)
  Quan <- reshape2::melt(Quan,c("Samples"))
  Quan$Samples <- factor( Quan$Samples, levels = unique(Quan$Samples))
  p3 <- ggplot(Quan,aes(x=variable, y = value, width=0.7)) +
    geom_bar(aes(fill=Samples),color="black",stat="identity", position=position_dodge()) +
    scale_fill_manual(values = ColorUse2) +
    # geom_text(aes(label = round(value,2)), size = 5) +
    theme_bw() +
    labs(x="", y="Percentage of each cell type in each sample") +
    theme(
      axis.title.x = element_text(color="black", size=16, face="bold"),
      axis.text.x = element_text(color="black", size=14, face="bold", angle = 60, hjust = 1),
      axis.title.y = element_text(color="black", size=16, face="bold"),
      axis.text.y = element_text(color="black", size=14, face="bold")
    )
  plot(p3)
  ggsave(paste0(OutPrefix,"cellNumbersPerCluster_dodged_",Sys.Date(),".pdf"), width = 8, height = 6)
}

##' Inference cell-cell communication on each sample using CellChat
##' @param Seurat_Obj A Seurat object
##' @param Seurat_MetaInfo Meta information of your Seurat object
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##' Run_CellChat(SeuratObj = iWAT_KO_HFD_data, Seurat_MetaInfo = iWAT_KO_HFD_meta, OutPrefix = "_iWAT_KO_HFD_CellChat_", ColorUse = iWAT_HFD_maps)
Run_CellChat <- function(SeuratObj, Seurat_MetaInfo, OutPrefix, ColorUse){
  
  if(file.exists(paste0(OutPrefix,".rds"))){
    cellchat <- readRDS(paste0(OutPrefix,".rds"))
  } else {
    cellchat <- createCellChat(object = SeuratObj, meta = Seurat_MetaInfo, group.by = "Anno")
    CellChatDB <- CellChatDB.mouse
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # dplyr::glimpse(CellChatDB$interaction)
    cellchat@DB <- CellChatDB.use
    # showDatabaseCategory(CellChatDB)
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 20) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
    cellchat <- projectData(cellchat, PPI.mouse)
    cellchat <- computeCommunProb(cellchat)  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05) #Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, paste0(OutPrefix,".rds"))
    future::plan("multisession", workers = 2) #Finish parallel
  }
  
  
  # cellchat <- updateCellChat(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))

  #Plot interactions and Weights
  pdf(paste0(OutPrefix, "_Interactions_",Sys.Date(),".pdf"), width = 16, height = 8)
  par(mfrow = c(1,3), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions", vertex.label.cex = 1.5,color.use = ColorUse)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 1.5,color.use = ColorUse)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,  title.name = "Interaction strength", vertex.label.cex = 1.5,color.use = ColorUse)
  dev.off()
  
  pdf(paste0(OutPrefix, "_Interactions_Seperate_",Sys.Date(),".pdf"), width = 12, height = 12)
  mat <- cellchat@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, label.edge = F, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = ColorUse)
  }
  dev.off()
  
  SigPathways <- cellchat@netP$pathways

  pdf(paste0(OutPrefix, "_PathwayInteractions_circle_",Sys.Date(),".pdf"), width = 12, height = 12)
  par(mfrow = c(4,4), xpd=TRUE)
  for (pathway in SigPathways){netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", color.use = ColorUse)}
  dev.off()

  pdf(paste0(OutPrefix, "_PathwayInteractions_chord_",Sys.Date(),".pdf"), width = 16, height = 16)
  par(mfrow = c(4,4), xpd=TRUE)
  for (pathway in SigPathways){netVisual_aggregate(cellchat, signaling = pathway, layout = "chord", color.use = ColorUse)}
  dev.off()
  
  pdf(paste0(OutPrefix, "_PathwayInteractions_netAnalysis_contribution_",Sys.Date(),".pdf"), width = 10, height = 10)
  for (pathway in SigPathways){p <- netAnalysis_contribution(cellchat, signaling = pathway);plot(p)}
  dev.off()
  
  # show all the significant interactions (L-R pairs)
  pdf(paste0(OutPrefix,"_PathwayInteractions_netVisual_bubble_",Sys.Date(),".pdf"), width = 12, height = 10)
  p <- netVisual_bubble(cellchat, remove.isolate = FALSE, color.text.use = ColorUse);plot(p)
  dev.off()
  
  pdf(paste0(OutPrefix,"_PathwayInteractions_PathwayGeneExpression_",Sys.Date(),".pdf"), width = 12, height = 12)
  for (pathway in SigPathways){p <- plotGeneExpression(cellchat, signaling = pathway, color.use = ColorUse, enriched.only = F);plot(p)}
  dev.off()

  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_network_",Sys.Date(),".pdf"), width = 8, height = 6)
  for (pathway in SigPathways){netAnalysis_signalingRole_network(cellchat, signaling = pathway, width = 8, height = 2.5, font.size = 10, color.use = ColorUse)}
  dev.off()

  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_scatter_",Sys.Date(),".pdf"), width = 6, height = 6)
  p <- netAnalysis_signalingRole_scatter(cellchat,color.use = ColorUse);plot(p)
  dev.off()
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  pdf(paste0(OutPrefix,"_PathwayInteractions_netAnalysis_signalingRole_heatmap_",Sys.Date(),".pdf"), width = 14, height = 10)
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = ColorUse, width = 10, height = 18, color.heatmap = "OrRd")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = ColorUse, width = 10, height = 18, color.heatmap = "OrRd")
  draw(ht1 + ht2)
  dev.off()

  pdf(paste0(OutPrefix, "_PathwayInteractions_Pattern_outgoing_",Sys.Date(),".pdf"), width = 8, height = 6)
  # selectK(cellchat, pattern = "outgoing")
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  p1 <- netAnalysis_river(cellchat, pattern = "outgoing",color.use = ColorUse)
  p2 <- netAnalysis_dot(cellchat, pattern = "outgoing",color.use = ColorUse)
  plot(p1+p2)
  dev.off()
  
  pdf(paste0(OutPrefix, "_PathwayInteractions_Pattern_incoming_",Sys.Date(),".pdf"), width = 8, height = 6)
  # selectK(cellchat, pattern = "incoming")
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  p1 <- netAnalysis_river(cellchat, pattern = "incoming",color.use = ColorUse)
  p2 <- netAnalysis_dot(cellchat, pattern = "incoming",color.use = ColorUse)
  plot(p1+p2)
  dev.off()
}

##' Comparison analysis of multiple datasets
##' @param CellChat_Obj_list A list of Seurat object
##' @param Seurat_MetaInfo Meta information of your Seurat object
##' @param ColorUse List of colors for plot
##' @param OutPrefix Prefix for output file
##' @examples
##' 
##' Run_CellChat_Compare(CellChat_Obj_list = iWAT_object.list, OutPrefix = "_iWAT_CellChat_Merged_", ColorUse = iWAT_HFD_maps)
Run_CellChat_Compare <- function(CellChat_Obj_list, OutPrefix, ColorUse){
  
  ListLen <- length(CellChat_Obj_list)
  combinations <- combn(ListLen,2)
  SigPathways <- c()
  for (i in 1:ListLen){SigPathways <- unique(c(SigPathways,CellChat_Obj_list[[i]]@netP$pathways))}
  
  cellchat <- mergeCellChat(CellChat_Obj_list, add.names = names(CellChat_Obj_list),cell.prefix = TRUE)
  # Compare the total number of interactions and interaction strength
  pdf(paste0(OutPrefix, "_compareInteractions_",Sys.Date(),".pdf"), width = 8, height = 8)
  gg1 <- compareInteractions(cellchat, show.legend = F)
  gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight")
  plot(gg1 + gg2)
  dev.off()

  # Compare the number of interactions and interaction strength among different cell populations
  pdf(paste0(OutPrefix, "_netVisual_diffInteraction_",Sys.Date(),".pdf"), width = 16, height = 10)
  for (i in 1:ncol(combinations)){
    SampleA <- combinations[1,i]
    SampleB <- combinations[2,i]
    par(mfrow = c(1,3), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T, color.use = ColorUse, comparison = c(SampleA, SampleB), label.edge = T,
                              title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
    netVisual_diffInteraction(cellchat, weight.scale = T, color.use = ColorUse, comparison = c(SampleA, SampleB), label.edge = F,
                              title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = ColorUse,comparison =  c(SampleA, SampleB))
  }
  dev.off()

  pdf(paste0(OutPrefix, "_netVisual_heatmap_",Sys.Date(),".pdf"), width = 12, height = 6)
  for (i in 1:ncol(combinations)){
    SampleA <- combinations[1,i]
    SampleB <- combinations[2,i]
    gg1 <- netVisual_heatmap(cellchat, color.use = ColorUse, comparison = c(SampleA, SampleB),
                             title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
    gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = ColorUse,comparison = c(SampleA, SampleB))
    plot(gg1 + gg2)
  }
  dev.off()

  pdf(paste0(OutPrefix, "_netVisual_circle_",Sys.Date(),".pdf"), width = 10, height = 10)
  weight.max <- getMaxWeight(CellChat_Obj_list, attribute = c("idents","count"))
  par(mfrow = c(3,2), xpd=TRUE)
  for (i in 1:length(CellChat_Obj_list)) {
    # netVisual_circle(CellChat_Obj_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], vertex.label.cex = 1.5,
    #                  edge.width.max = 12, title.name = paste0("Number of interactions - ", names(CellChat_Obj_list)[i]), color.use = ColorUse)
    netVisual_circle(CellChat_Obj_list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2],vertex.label.cex = 1.5,edge.label.cex = 1.2,
                     edge.width.max = 12, title.name = paste0("Number of interactions - ", names(CellChat_Obj_list)[i]), color.use = ColorUse)
  }
  dev.off()

  # Identify and visualize the conserved and context-specific signaling pathways
  # Compare the overall information flow of each signaling pathway
  pdf(paste0(OutPrefix, "_rankNet_",Sys.Date(),".pdf"), width = 12, height = 6)
  for (i in 1:ncol(combinations)){
    SampleA <- combinations[1,i]
    SampleB <- combinations[2,i]
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(SampleA, SampleB))
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(SampleA, SampleB))
    plot(gg1 + gg2)
  }
  dev.off()

  # Compare outgoing (or incoming) signaling associated with each cell population
  pdf(paste0(OutPrefix,"_netAnalysis_signalingRole_heatmap_",Sys.Date(),".pdf"), width = 14, height = 12)
  for (i in 1:ListLen){
    par(mfrow = c(2,2), xpd=TRUE)
    ht <- netAnalysis_signalingRole_heatmap(CellChat_Obj_list[[i]], pattern = "all", signaling = SigPathways, color.use = ColorUse,
                                            title = names(CellChat_Obj_list)[i], width = 10, height = 18, color.heatmap = "OrRd")
    draw(ht)
  }
  dev.off()
  
  # Identify the upgulated and down-regulated signaling ligand-receptor pairs
  pdf(paste0(OutPrefix, "_netVisual_bubble_",Sys.Date(),".pdf"), width = 14, height = 12)
  for (i in 1:ncol(combinations)){
    SampleA <- combinations[1,i]
    SampleB <- combinations[2,i]
    SigPathways <- intersect(CellChat_Obj_list[[SampleA]]@netP$pathways,CellChat_Obj_list[[SampleB]]@netP$pathways)
    p <- netVisual_bubble(cellchat, comparison = c(SampleA, SampleB), signaling = SigPathways, angle.x = 45,
                          title.name = paste0(names(CellChat_Obj_list[SampleA]),"_vs_",names(CellChat_Obj_list[SampleB])))
    plot(p)
  }
  dev.off()
  
  pdf(paste0(OutPrefix, "_plotGeneExpression_",Sys.Date(),".pdf"), width = 12, height = 12)
  cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = unique(cellchat@meta$orig.ident)) # set factor level
  for (pathway in SigPathways){
    p <- plotGeneExpression(cellchat, signaling = pathway, split.by = "datasets", colors.ggplot = T,color.use = ColorUse)
    plot(p)}
  dev.off()
}

