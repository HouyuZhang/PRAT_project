source("0.ScFunctions.R")

#===============================================================================
# 1. mPRAT-iWAT ACE integration
#===============================================================================
mPRATiWAT_prefix <- "mPRATiWAT_CE_"
if(file.exists(paste0(mPRATiWAT_prefix,"Seurat.rds"))){
  mPRATiWAT_combined <- readRDS(paste0(mPRATiWAT_prefix,"Seurat.rds"))
} else {
  mPRAT_ACE_data <- ReadCB_h5("../_rawdata/mPRAT-6mo-ACE/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iWAT_ACE_data <- ReadCB_h5("../_rawdata/iWAT-6mo-ACE/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  mPRAT_ACE_Obj <- CreateSeuratObject(counts = mPRAT_ACE_data, project = "mPRAT-ACE", min.cells = 10)
  iWAT_ACE_Obj <- CreateSeuratObject(counts = iWAT_ACE_data, project = "iWAT-ACE", min.cells = 10)
  mPRATiWAT_combined <- RunSeurat(ObjList = c(mPRAT_ACE_Obj,iWAT_ACE_Obj), Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), 
                                 OutPrefix = mPRATiWAT_prefix, NormalizationMethod = "VST", MitoPercent = 15)
  
  Idents(mPRATiWAT_combined) <- mPRATiWAT_combined$integrated_snn_res.0.1
  DimPlot(mPRATiWAT_combined, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = mPRATiWAT_combined, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRATiWAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWAT_combined, ColorUse = hughie_color, OutPrefix = mPRATiWAT_prefix)
}

mPRATiWAT_combined <- RenameIdents(mPRATiWAT_combined,  `0` = "Adipocyte", `1` = "Adipocyte", `2` = "Adipocyte",`3` = "ASPC", `4` = "Macrophage", 
                                   `5` = "Endothelial Cell",`6` = "Mesothelial Cell",`7` = "Pericyte", `8` = "Macrophage")
mPRATiWAT_combined$Myannotation <- mPRATiWAT_combined@active.ident
mPRATiWAT_cellRanks <- c("ASPC","Adipocyte","Macrophage","Endothelial Cell","Pericyte","Mesothelial Cell")
mPRATiWAT_color <- c("#0072B5FF","#BC3C29FF", "#7E6148FF","#E18727FF","#EE4C97FF", "#631879FF")
names(mPRATiWAT_color) <- mPRATiWAT_cellRanks
Idents(mPRATiWAT_combined) <- factor(Idents(mPRATiWAT_combined), levels= mPRATiWAT_cellRanks)
mPRATiWAT_combined$orig.ident <- factor(mPRATiWAT_combined$orig.ident, levels = c("mPRAT-ACE","iWAT-ACE"))
mPRATiWAT_CE_SAMcolors <- c("black","#EE4C97FF")

PlotGenes <- c("Pdgfra","Nova1", #ASPC
               "Plin1","Adipoq","Acvr1c",#Adipocytes
               "Ptprc","Mrc1", #Macrophage
               "Cyyr1","Flt1", #Endothelial
               "Abcc9","Rgs5", #Pericyte
               "Muc16", "Msln" #Mesothelial
)

(p1 <- VlnPlot(mPRATiWAT_combined, features = PlotGenes, fill.by = "ident", cols = mPRATiWAT_color, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRATiWAT_combined, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRATiWAT_SAMcolors, assay = "RNA",split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRATiWAT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width= 14, height = 8)

(p1 <- DimPlot(mPRATiWAT_combined, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mPRATiWAT_color, label.color = "white"))
(p2 <- DimPlot(mPRATiWAT_combined, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRATiWAT_SAMcolors))
p1 + p2
ggsave(paste0(mPRATiWAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mPRATiWAT_combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRATiWAT_color, label.color = "white", order = T)
ggsave(paste0(mPRATiWAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = subset(mPRATiWAT_combined, idents = mPRATiWAT_cellRanks), OutPrefix = mPRATiWAT_prefix, ColorUse1 = mPRATiWAT_color, ColorUse2 = mPRATiWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWAT_combined, ColorUse = mPRATiWAT_color, OutPrefix = mPRATiWAT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRATiWAT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = mPRATiWAT_cellRanks)

#===============================================================================
# 2. mPRAT-iWAT ACE Adipocytes
#===============================================================================
mPRATiWAT_CE_Ad_prefix <- "mPRATiWAT_CE_Adipocyte_"
if(file.exists(paste0(mPRATiWAT_CE_Ad_prefix, "Seurat.rds"))){
  mPRATiWAT_CE_Ad_comb <- readRDS(paste0(mPRATiWAT_CE_Ad_prefix, "Seurat.rds"))
} else {
  mPRATiWAT_CE_Ad_merged <- ScaleData(subset(mPRATiWAT_combined, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  clustree(mPRATiWAT_CE_Ad_merged)
  ggsave(paste0(mPRATiWAT_CE_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(mPRATiWAT_CE_Ad_merged) <- mPRATiWAT_CE_Ad_merged$integrated_snn_res.0.2
  DimPlot(mPRATiWAT_CE_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRATiWAT_CE_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRATiWAT_CE_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWAT_CE_Ad_merged, ColorUse = hughie_color, OutPrefix = mPRATiWAT_CE_Ad_prefix)
}

mPRATiWAT_cells <- c("mPRAT-ad1", "mPRAT-ad3","mPRAT-ad4", "iWAT-ad1", "iWAT-ad2", "iWAT-ad3","iWAT-ad4","iWAT-ad5")
mPRATiWAT_color <- c("#BC3C29FF", "#B2DF8A","#20854EFF","#BC3C29FF","#FF7F00","#0072B5FF","#6A3D9A","#20854EFF")
names(mPRATiWAT_cells) <- mPRATiWAT_color

mPRATiWAT_CE_Ad_comb$AllID <- mPRATiWAT_CE_Ad_comb$mPRATCE_anno
mPRATiWAT_CE_Ad_comb@meta.data[mPRATiWAT_CE_Ad_comb$orig.ident == "iWAT-ACE",]$AllID <- mPRATiWAT_CE_Ad_comb$iWATCE_anno[mPRATiWAT_CE_Ad_comb$orig.ident == "iWAT-ACE"]
mPRATiWAT_CE_Ad_comb$AllID <- factor(mPRATiWAT_CE_Ad_comb$AllID, levels = mPRATiWAT_cells)
Idents(mPRATiWAT_CE_Ad_comb) <- mPRATiWAT_CE_Ad_comb$integrated_snn_res.0.3

(p1 <- DimPlot(mPRATiWAT_CE_Ad_comb, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = mPRATiWAT_CE_SAMcolors))
(p2 <- DimPlot(mPRATiWAT_CE_Ad_comb, reduction = "umap", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = hughie_color))
(p3 <- DimPlot(mPRATiWAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", group.by = "AllID",ncol = 2, label = F, repel = T, label.box =T, order = T, cols = mPRATiWAT_color))
(p1|p2)/p3
ggsave(paste0(mPRATiWAT_CE_Ad_prefix,"Compare_CE",Sys.Date(),".pdf"), width = 13, height = 13)

Idents(mPRATiWAT_CE_Ad_comb) <- mPRATiWAT_CE_Ad_comb$AllID
PlotGenes <- c("Acvr1c","Ucp1","Cidea","Cyp2e1","Slit3","Slc7a10","Aldh1a1","Lep")
VlnPlot(mPRATiWAT_CE_Ad_comb, features = PlotGenes, fill.by = "ident", cols = mPRATiWAT_color, stack = T, flip=T, assay = "RNA")+ NoLegend()
ggsave(paste0(mPRATiWAT_CE_Ad_prefix,"Compare_CE_VlnPlot_", Sys.Date(), ".pdf"), width = 6, height = 8)

TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWAT_CE_Ad_comb, ColorUse = mPRATiWAT_color, OutPrefix = mPRATiWAT_CE_Ad_prefix, NMarkers = 20)
expression_averge <- AverageExpression(mPRATiWAT_CE_Ad_comb)

markers <- read.csv("mPRATiWAT_CE_Adipocyte_Allmarkers.csv") %>% filter(p_val_adj < 0.01)
pdf(paste0(mPRATiWAT_CE_Ad_prefix,"compare_adipocytes_DEG_heatmap_", Sys.Date(),".pdf"), width = 8, height = 12)
pheatmap::pheatmap(expression_averge$RNA[unique(markers$gene),], show_colnames = T, show_rownames = T, scale = "row")
dev.off()

#------------------- label transfer from Time-series dataset ---------------------------
mi_Transfer_prefix <- "mPRAT_iWAT_CE_Adipocyte_Transfer_"
if(file.exists(paste0(mi_Transfer_prefix,"Seurat.rds"))){
  mPRATiWAT_CE_Ad_comb <- readRDS(paste0(mi_Transfer_prefix,"Seurat.rds"))
} else {
  mPRATiWAT_CE_Ad_comb <- subset(mPRATiWAT_CE_Ad_comb, subset = orig.ident == "mPRAT-6mo-ACE")
  mPRAT_Ad_merged <- readRDS(paste0("../1. mPRAT_126mo/mPRAT_Adipocytes_Seurat.rds"))
  mPRAT_Ad_merged6mo <- subset(mPRAT_Ad_merged, subset = orig.ident == "mPRAT-6mo")

  mPRATanchors <- FindTransferAnchors(reference = mPRAT_Ad_merged, query = mPRATiWAT_CE_Ad_comb, dims = 1:15, reference.reduction = "pca")
  predictions <- TransferData(anchorset = mPRATanchors, refdata = mPRAT_Ad_merged$Myannotation_Ad, dims = 1:15)
  mPRATiWAT_CE_Ad_comb <- AddMetaData(mPRATiWAT_CE_Ad_comb, metadata = predictions)
  mPRATiWAT_CE_Ad_comb$predicted.id <- factor(mPRATiWAT_CE_Ad_comb$predicted.id, levels = c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4"))
  saveRDS(mPRATiWAT_CE_Ad_comb, file = paste0(Transfer_prefix,"Seurat.rds"))
  
  mPRAT_CE_Ad_celltypeRanks <- c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4")
  mPRAT_CE_Ad_colors <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")
  names(mPRAT_CE_Ad_colors) <- mPRAT_CE_Ad_celltypeRanks
  (p <- DimPlot(mPRATiWAT_CE_Ad_comb, reduction = "umap", group.by = "predicted.id", split.by = "orig.ident", label = TRUE, label.size = 3, 
                repel = TRUE, cols = mPRAT_CE_Ad_colors) + NoLegend())
  ggsave(paste0(mi_Transfer_prefix, "UMAP_2", Sys.Date(), ".pdf"), width = 8, height = 10)
}
