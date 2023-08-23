source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================
PlotObjMetrices(Seurat_Obj = iWAT_merged, OutPrefix = iWAT_prefix, ColorUse = iWAT_SAMcolors)
GetCellComposition(iWAT_merged)

#===============================================================================
# 1. iWAT RT CE integration
#===============================================================================
iWAT_prefix <- "iWAT_RTvsACE_"
if(file.exists(paste0(iWAT_prefix,"Seurat.rds"))){
  iWAT_merged <- readRDS(paste0(iWAT_prefix,"Seurat.rds"))
} else {
  iWAT_data <- ReadCB_h5("../_rawdata/iWAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iWAT_ACE_data <- ReadCB_h5("../_rawdata/iWAT-6mo-ACE/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iWAT_Obj <- CreateSeuratObject(counts = iWAT_data, project = "iWAT-RT", min.cells = 10)
  iWAT_ACE_Obj <- CreateSeuratObject(counts = iWAT_ACE_data, project = "iWAT-ACE", min.cells = 10)
  iWAT_merged <- RunSeurat(ObjList = c( iWAT_Obj, iWAT_ACE_Obj),Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), 
                                 OutPrefix = iWAT_prefix, NormalizationMethod = "VST", MitoPercent = 15)
  
  Idents(iWAT_merged) <- iWAT_merged$integrated_snn_res.0.1
  DimPlot(iWAT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = iWAT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iWAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_merged, ColorUse = hughie_color, OutPrefix = iWAT_prefix)
}

iWAT_merged <- RenameIdents(iWAT_merged,  `0` = "Adipocyte", `1` = "ASPC", `2` = "Adipocyte",`3` = "Macrophage", `4` = "T Cell", 
                                  `5` = "Endothelial Cell", `6` = "B Cell",`7` = "Adipocyte", `8` = "Dendritic Cell", `9` = "Pericyte", `10` = "Neutrophill")

iWAT_merged$Myannotation <- iWAT_merged@active.ident
iWAT_cellRanks <- c("ASPC","Adipocyte", "Macrophage","Neutrophill","Dendritic Cell","T Cell","B Cell","Endothelial Cell","Pericyte")
iWAT_color <- c("#0072B5FF","#BC3C29FF", "#7E6148FF", "#B09C85FF", "#808180FF","#E18727FF","#EE4C97FF","#631879FF","#3B4992FF")
names(iWAT_color) <- iWAT_cellRanks
Idents(iWAT_merged) <- factor(Idents(iWAT_merged), levels= iWAT_cellRanks)
iWAT_merged$orig.ident <- factor(iWAT_merged$orig.ident, levels = c("iWAT-RT","iWAT-ACE"))
iWAT_SAMcolors <- c("#BC3C29FF","#0072B5FF")

PlotGenes <- c("Pdgfra","Nova1", #ASPC
               "Plin1","Acvr1c",#Adipocytes
               "Dock2","Ptprc","Adgre1","Mrc1", #Macrophage
               "Csf3r", #Neutrophill
               "Flt3", #Dendritic Cell
               "Prkcq","Themis",# T Cell
               "Pax5","Bcl11a", # B Cell
               "Cyyr1","Flt1", #Endothelial
               "Abcc9","Rgs5"#Pericyte
)

(p1 <- VlnPlot(iWAT_merged, features = PlotGenes, fill.by = "ident", cols = iWAT_color, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iWAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iWAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iWAT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width= 14, height = 8)

(p1 <- DimPlot(iWAT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = iWAT_color, label.color = "white"))
(p2 <- DimPlot(iWAT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iWAT_SAMcolors))
p1 + p2
ggsave(paste0(iWAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(iWAT_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iWAT_color, label.color = "white", order = T)
ggsave(paste0(iWAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 14, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iWAT_merged, OutPrefix = iWAT_prefix, ColorUse1 = iWAT_color, ColorUse2 = hughie_color)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_merged, ColorUse = iWAT_color, OutPrefix = iWAT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iWAT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = iWAT_cellRanks)
saveRDS(iWAT_merged, file = paste0(iWAT_prefix,"Seurat.rds"))

#===============================================================================
# 2. iWAT RT CE Adipocytes
#===============================================================================
iWAT_Ad_prefix <- "iWAT_RTvsACE_Adipocytes_"
if(file.exists(paste0(iWAT_Ad_prefix,"Seurat.rds"))){
  iWAT_Ad_merged <- readRDS(paste0(iWAT_Ad_prefix,"Seurat.rds"))
} else {
  iWAT_Ad_merged <- ScaleData(subset(iWAT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(iWAT_Ad_merged)
  ggsave(paste0(iWAT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)

  Idents(iWAT_Ad_merged) <- iWAT_Ad_merged$integrated_snn_res.0.15
  DimPlot(iWAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iWAT_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iWAT_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_Ad_merged, ColorUse = hughie_color, OutPrefix = iWAT_Ad_prefix)
}

iWAT_Ad_merged <- RenameIdents(iWAT_Ad_merged, `0` = "iWAT-ad4", `1` = "iWAT-ad5", `2` = "iWAT-ad2",`3` = "iWAT-ad3",`4` = "iWAT-ad1")
iWAT_Ad_merged$Myannotation_Ad <- iWAT_Ad_merged@active.ident
iWAT_AdRanks <- c("iWAT-ad1","iWAT-ad2","iWAT-ad3","iWAT-ad4","iWAT-ad5")
iWAT_Ad_colors <- c("#BC3C29FF","#FF7F00","#0072B5FF","#6A3D9A","#20854EFF")
names(iWAT_Ad_colors) <- iWAT_AdRanks
Idents(iWAT_Ad_merged) <- factor(Idents(iWAT_Ad_merged), levels = iWAT_AdRanks)

PlotGenes <- c("Plin1","Ucp1","Cidea","Acss2","Fasn","Acly",
               "Lep","Cyp2e1","Aldh1a1","Slc7a10","Slit3","Npr3"
)

(p1 <- VlnPlot(iWAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = iWAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iWAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iWAT_SAMcolors, assay = "RNA", split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iWAT_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 16, height = 12)

(p1 <- DimPlot(iWAT_Ad_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = iWAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(iWAT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iWAT_SAMcolors))
p1 + p2
ggsave(paste0(iWAT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 12, height = 5)

DimPlot(iWAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iWAT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(iWAT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iWAT_Ad_merged, OutPrefix = iWAT_Ad_prefix, ColorUse1 = iWAT_Ad_colors, ColorUse2 = iWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_Ad_merged, ColorUse = iWAT_Ad_colors, OutPrefix = iWAT_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iWAT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = iWAT_AdRanks)

FeaturePlot(iWAT_Ad_merged, features =  c("Ucp1","Acly","Cyp2e1","Lep"), min.cutoff = "q10", order = F, pt.size = 0.3)
ggsave(paste0(iWAT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(iWAT_Ad_merged, file = paste0(iWAT_Ad_prefix,"Seurat.rds"))

#===============================================================================
# 2. iWAT RT CE ASPC
#===============================================================================
iWAT_ASPC_prefix <- "iWAT_RTvsACE_ASPC_"
if(file.exists(paste0(iWAT_ASPC_prefix,"Seurat.rds"))){
  iWAT_ASPC_merged <- readRDS(paste0(iWAT_ASPC_prefix,"Seurat.rds"))
} else {
  iWAT_ASPC_merged <- ScaleData(subset(iWAT_merged, idents = c("ASPC"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  clustree(iWAT_ASPC_merged)
  ggsave(paste0(iWAT_ASPC_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)

  Idents(iWAT_ASPC_merged) <- iWAT_ASPC_merged$integrated_snn_res.0.2
  DimPlot(iWAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iWAT_ASPC_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iWAT_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_ASPC_merged, ColorUse = hughie_color, OutPrefix = iWAT_ASPC_prefix)
}

iWAT_ASPC_merged <- RenameIdents(iWAT_ASPC_merged, `0` = "iWAT-aspc2", `1` = "iWAT-aspc1", `2` = "iWAT-aspc3")
iWAT_ASPC_merged$Myannotation_ASPC <- iWAT_ASPC_merged@active.ident
iWAT_ASPCRanks <- c("iWAT-aspc1","iWAT-aspc2","iWAT-aspc3")
iWAT_ASPC_colors <- c("#0072B5FF","#E18727FF", "#6A3D9A")
names(iWAT_ASPC_colors) <- iWAT_ASPCRanks
Idents(iWAT_ASPC_merged) <- factor(Idents(iWAT_ASPC_merged), levels = iWAT_ASPCRanks)

PlotGenes <- c("Pdgfra","Dpp4","Anxa3","Pdgfrb","Pparg","Cd36","Fabp4")
(p1 <- VlnPlot(iWAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = iWAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iWAT_ASPC_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iWAT_SAMcolors, assay = "RNA", split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iWAT_ASPC_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 16, height = 12)

(p1 <- DimPlot(iWAT_ASPC_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = iWAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(iWAT_ASPC_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iWAT_SAMcolors))
p1 + p2
ggsave(paste0(iWAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 16, height = 8)

DimPlot(iWAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iWAT_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(iWAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = iWAT_ASPC_merged, OutPrefix = iWAT_ASPC_prefix, ColorUse1 = iWAT_ASPC_colors, ColorUse2 = iWAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iWAT_ASPC_merged, ColorUse = iWAT_ASPC_colors, OutPrefix = iWAT_ASPC_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iWAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4,cellRanks = iWAT_ASPCRanks)

FeaturePlot(iWAT_ASPC_merged, features =  c("Pdgfra","Dpp4","Pdgfrb","Fabp4"), min.cutoff = "q9", order = T, pt.size=0.3)
ggsave(paste0(iWAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 8, height = 10)
saveRDS(iWAT_ASPC_merged, file = paste0(iWAT_ASPC_prefix,"Seurat.rds"))
