source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================
PlotObjMetrices(Seurat_Obj = subset(iBAT_CE_comb, subset = orig.ident == "iBAT-6mo-ACE"), OutPrefix = iBAT_CE_prefix, ColorUse = iBAT_SAMcolors[2])
GetCellComposition(iBAT_CE_comb)

#===============================================================================
# 1. iBAT RT+ACE
#===============================================================================
iBAT_CE_prefix <- "iBAT_CE_ACells_"
if(file.exists(paste0(iBAT_CE_prefix,"Seurat.rds"))){
  iBAT_CE_comb <- readRDS(paste0(iBAT_CE_prefix,"Seurat.rds"))
} else {
  iBAT6mo_data <-    ReadCB_h5("../_rawdata/iBAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = T)
  iBAT6mo_CE_data <- ReadCB_h5("../_rawdata/iBAT-6mo-ACE/CellBender_feature_bc_matrix_filtered.h5", use.names = T)
  iBAT6mo_Obj <- CreateSeuratObject(counts = iBAT6mo_data, project = "iBAT-6mo-RT", min.cells = 10)
  iBAT6mo_CE_Obj <- CreateSeuratObject(counts = iBAT6mo_CE_data, project = "iBAT-6mo-ACE", min.cells = 10)
  iBAT_CE_comb <- RunSeurat(ObjList = c(iBAT6mo_Obj, iBAT6mo_CE_Obj), OutPrefix = iBAT_CE_prefix,
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)
  
  Idents(iBAT_CE_comb) <- iBAT_CE_comb$integrated_snn_res.0.1
  DimPlot(iBAT_CE_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = T, repel = T, label.box =F, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iBAT_CE_comb, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBAT_CE_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_CE_comb, ColorUse = hughie_color, OutPrefix = iBAT_CE_prefix)
}

iBAT_CE_comb <- RenameIdents(iBAT_CE_comb,  `0` = "Adipocyte", `1` = "Adipocyte", `2` = "Adipocyte",`3` = "ASPC", `4` = "Endothelial Cell", 
                             `5` = "Macrophage", `6` = "Adipocyte",  `7` = "Pericyte")
iBAT_CE_celltypeRanks <- c("ASPC","Adipocyte","Endothelial Cell","Pericyte","Macrophage")
iBAT_CE_colorMaps <- c("#0072B5FF","#BC3C29FF","#631879FF","#3B4992FF","#7E6148FF")
names(iBAT_CE_colorMaps) <- iBAT_CE_celltypeRanks
iBAT_CE_comb$Myannotation <- iBAT_CE_comb@active.ident
Idents(iBAT_CE_comb) <- factor(Idents(iBAT_CE_comb), levels = iBAT_CE_celltypeRanks)
iBAT_CE_comb$orig.ident <- factor(iBAT_CE_comb$orig.ident, levels = c("iBAT-6mo-RT", "iBAT-6mo-ACE"))
iBAT_CE_SAMcolors <- c("#BC3C29FF","#0072B5FF")

PlotGenes <- c("Pdgfra","Nova1", # ASPC
               "Plin1","Ucp1", # Adipocyte
               "Cyyr1","Flt1", # Endothelial Cell
               "Rgs5","Abcc9", # Pericyte
               "Dock2","Mrc1"# Macrophage
)

(p1 <- VlnPlot(iBAT_CE_comb, features = PlotGenes, fill.by = "ident", cols = iBAT_CE_colorMaps, stack = T, flip=T,assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iBAT_CE_comb, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = hughie_color, assay = "RNA")) + FontSize(x.text = 16, y.text = 16)
p1 + p2
ggsave(paste0(iBAT_CE_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 20, height = 12)

(p1 <- DimPlot(iBAT_CE_comb, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = iBAT_CE_colorMaps, label.color = "white"))
(p2 <- DimPlot(iBAT_CE_comb, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = hughie_color))
p1 + p2
ggsave(paste0(iBAT_CE_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 16, height = 8)

DimPlot(iBAT_CE_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iBAT_CE_colorMaps, label.color = "white", order = T)
ggsave(paste0(iBAT_CE_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = subset(iBAT_CE_comb, ident = iBAT_CE_celltypeRanks), OutPrefix = iBAT_CE_prefix, ColorUse1 = iBAT_CE_colorMaps, ColorUse2 = iBAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_CE_comb, ColorUse = iBAT_CE_colorMaps, OutPrefix = iBAT_CE_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iBAT_CE_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = iBAT_CE_celltypeRanks)

FeaturePlot(iBAT_CE_comb, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), split.by = "orig.ident", min.cutoff = "q10", order = T, pt.size=0.3) 
ggsave(paste0(iBAT_CE_prefix,"FeaturePlot_", Sys.Date(), ".pdf"), width = 10, height = 18)
saveRDS(iBAT_CE_comb, file = paste0(iBAT_CE_prefix,"Seurat.rds"))

#===============================================================================
# 2. iBAT RT+ACE Adipocytes
#===============================================================================
iBATCE_Ad_prefix <- "iBAT_CE_Adipocyte_"
if(file.exists(paste0(iBATCE_Ad_prefix,"Seurat.rds"))){
  iBAT_CE_Ad_comb <- readRDS(paste0(iBATCE_Ad_prefix,"Seurat.rds"))
} else {
  iBAT_CE_Ad_comb <- subset(iBAT_CE_comb, idents = c("Adipocyte")) %>% ScaleData() %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(iBAT_CE_Ad_comb)
  ggsave(paste0(iBATCE_Ad_prefix, "clustree_", Sys.Date(), ".pdf"), width = 12, height = 10)
  saveRDS(iBAT_CE_Ad_comb, file = paste0(iBATCE_Ad_prefix,"Seurat.rds"))
  
  Idents(iBAT_CE_Ad_comb) <- iBAT_CE_Ad_comb$integrated_snn_res.0.2
  DimPlot(iBAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iBAT_CE_Ad_comb, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBATCE_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_CE_Ad_comb, ColorUse = hughie_color, OutPrefix = iBATCE_Ad_prefix, NMarkers = 20)
  
  iBAT_CE_Ad_comb <- RenameIdents(iBAT_CE_Ad_comb, `0` = "iBAT-ad1", `1` = "iBAT-ad3", `2` = "iBAT-ad1",`3` = "iBAT-ad1",`4` = "iBAT-ad2",`5` = "iBAT-ad4")
  iBAT_Ad_celltypeRanks <- c("iBAT-ad1","iBAT-ad2","iBAT-ad3","iBAT-ad4")
  iBAT_Ad_colors <- c("#BC3C29FF", "#EE4C97FF","#FB9A99","#20854EFF")
  names(iBAT_Ad_colors) <- iBAT_Ad_celltypeRanks
  Idents(iBAT_CE_Ad_comb) <- factor(Idents(iBAT_CE_Ad_comb), levels = iBAT_Ad_celltypeRanks)
  iBAT_CE_Ad_comb$ACEAnno <- iBAT_CE_Ad_comb@active.ident
  iBAT_CE_Ad_comb$orig.ident <- factor(iBAT_CE_Ad_comb$orig.ident, levels = c("iBAT-6mo-RT", "iBAT-6mo-ACE"))
}

PlotGenes <- c("Plin1","Ucp1","Gk","Cidea","Cyp2e1","Aldh1a1","Slc7a10")
(p1 <- VlnPlot(iBAT_CE_Ad_comb, features = PlotGenes, fill.by = "ident", cols = iBAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iBAT_CE_Ad_comb, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iBAT_CE_SAMcolors, assay = "RNA", split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iBATCE_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width =10, height = 8)

(p1 <- DimPlot(iBAT_CE_Ad_comb, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T,  cols = iBAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(iBAT_CE_Ad_comb, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iBAT_CE_SAMcolors))
p1 + p2
ggsave(paste0(iBATCE_Ad_prefix, "UMAP-ByClusters_", Sys.Date(), ".pdf"), width = 14, height = 6)

DimPlot(iBAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iBAT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(iBATCE_Ad_prefix,"UMAP-Byorig.ident_", Sys.Date(), ".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iBAT_CE_Ad_comb, OutPrefix = iBATCE_Ad_prefix, ColorUse1 = iBAT_Ad_colors, ColorUse2 = iBAT_CE_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_CE_Ad_comb, ColorUse = iBAT_Ad_colors, OutPrefix = iBATCE_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iBATCE_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = iBAT_Ad_celltypeRanks)

FeaturePlot(iBAT_CE_Ad_comb, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), ncol = 2, split.by = "orig.ident", min.cutoff = "q10", order = T, pt.size=0.3) 
ggsave(paste0(iBATCE_Ad_prefix,"FeaturePlot_", Sys.Date(), ".pdf"), width = 10, height = 18)

#------------------- Compare annotation with time-series dataset ---------------------------
iBAT_Ad_combined1 <- readRDS(paste0("../2. iBAT_126mo/iBAT_Adipocytes_Seurat.rds"))
TS_metadata <- subset(iBAT_Ad_combined1, subset = orig.ident == "iBAT-6mo")@meta.data
rownames(TS_metadata) <- gsub("_3","",rownames(TS_metadata))
TS_metadata$Myannotation_Ad <- factor(TS_metadata$Myannotation_Ad,levels = c("iBAT-ad1","iBAT-ad2","iBAT-ad3","iBAT-ad4"))

iBAT_CE_Ad_combRT <- subset(iBAT_CE_Ad_comb, subset = orig.ident == "iBAT-6mo-RT")
Cells <- gsub("_1","",rownames(iBAT_CE_Ad_combRT@meta.data))
iBAT_CE_Ad_combRT$TS_anno <- TS_metadata[Cells,]$Myannotation_Ad

Idents(iBAT_CE_Ad_combRT) <- iBAT_CE_Ad_combRT$TS_anno
p1 <- DimPlot(iBAT_CE_Ad_combRT, reduction = "umap",ncol = 2, label = F, repel = T, label.box =T, order = T, cols = iBAT_CE_Ad_colors)
Idents(iBAT_CE_Ad_combRT) <- iBAT_CE_Ad_combRT$predicted.id
p2 <- DimPlot(iBAT_CE_Ad_comb, reduction = "umap",split.by = "orig.ident",ncol = 2, label = F, repel = T, label.box =T, order = T, cols = iBAT_CE_Ad_colors)
(p1)/p2
ggsave(paste0(iBATCE_Ad_prefix,"Compare_RTACE",Sys.Date(),".pdf"), width = 10, height = 10)

iBAT_CE_Ad_combRT@meta.data[,c("TS_anno","predicted.id")] %>% dplyr::count(TS_anno, predicted.id) %>% filter(TS_anno != "NA") %>% 
  ggplot(aes(axis1 = TS_anno, axis2 = predicted.id, y = n)) +
  geom_alluvium(aes(fill = TS_anno),alpha=0.5) +
  geom_stratum(aes(fill = TS_anno)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("#BC3C29FF", "#EE4C97FF","#FB9A99","#33A02C")) + 
  theme_void() + theme(legend.position = "none")
ggsave(paste0(iBATCE_Ad_prefix,"SankeyPlot_",Sys.Date(),".pdf"), width = 4, height = 8)

#===============================================================================
# 3. iBAT RT+ACE ASPC
#===============================================================================
iBATCE_ASPC_prefix <- "iBAT_CE_ASPC_"
if(file.exists(paste0(iBATCE_ASPC_prefix,"Seurat.rds"))){
  iBATCE_ASPC_comb <- readRDS(paste0(iBATCE_ASPC_prefix,"Seurat.rds"))
} else {
  iBATCE_ASPC_comb <- subset(iBAT_CE_comb, idents = c("ASPC")) %>% ScaleData() %>% RunPCA(npcs = 30) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>% FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(iBATCE_ASPC_comb)
  ggsave(paste0(iBATCE_ASPC_prefix, "clustree_", Sys.Date(), ".pdf"), width = 12, height = 10)
  
  Idents(iBATCE_ASPC_comb) <- iBATCE_ASPC_comb$integrated_snn_res.0.1
  DimPlot(iBATCE_ASPC_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iBATCE_ASPC_comb, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBATCE_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBATCE_ASPC_comb, ColorUse = hughie_color, OutPrefix = iBATCE_ASPC_prefix,NMarkers = 20)
}

iBATCE_ASPC_comb <- RenameIdents(iBATCE_ASPC_comb, `0` = "iBAT-aspc2", `1` = "iBAT-aspc3", `2` = "iBAT-aspc1")
iBAT_ASPC_celltypeRanks <- c("iBAT-aspc1","iBAT-aspc2","iBAT-aspc3")
iBAT_ASPC_colors <- c("#1F78B4","#FDBF6F","#CAB2D6")
names(iBAT_ASPC_colors) <- iBAT_ASPC_celltypeRanks
Idents(iBATCE_ASPC_comb) <- factor(Idents(iBATCE_ASPC_comb), levels = iBAT_ASPC_celltypeRanks)
iBATCE_ASPC_comb$ACEAnno <- iBATCE_ASPC_comb@active.ident
iBATCE_ASPC_comb$orig.ident <- factor(iBATCE_ASPC_comb$orig.ident, levels = c("iBAT-6mo-RT", "iBAT-6mo-ACE"))

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36","Cidea")
VlnPlot(iBATCE_ASPC_comb, features = PlotGenes, fill.by = "ident", cols = hughie_color, stack = T, flip=T,assay = "RNA")

(p1 <- VlnPlot(iBATCE_ASPC_comb, features = PlotGenes, fill.by = "ident", cols = iBAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iBATCE_ASPC_comb, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iBAT_SAMcolors, assay = "RNA",split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iBATCE_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width =10, height = 8)

(p1 <- DimPlot(iBATCE_ASPC_comb, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T,  cols = iBAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(iBATCE_ASPC_comb, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iBAT_SAMcolors))
p1 + p2
ggsave(paste0(iBATCE_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(), ".pdf"), width = 14, height = 6)

DimPlot(iBATCE_ASPC_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = iBAT_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(iBATCE_ASPC_prefix,"UMAP-Byorig.ident_", Sys.Date(), ".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iBATCE_ASPC_comb, OutPrefix = iBATCE_ASPC_prefix, ColorUse1 = iBAT_ASPC_colors, ColorUse2 = iBAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBATCE_ASPC_comb, ColorUse = iBAT_ASPC_colors, OutPrefix = iBATCE_ASPC_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iBATCE_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = iBAT_ASPC_celltypeRanks)

FeaturePlot(iBATCE_ASPC_comb, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), ncol = 2, split.by = "orig.ident", min.cutoff = "q10", order = T, pt.size=0.3) 
ggsave(paste0(iBATCE_ASPC_prefix,"FeaturePlot_", Sys.Date(), ".pdf"), width = 10, height = 18)

saveRDS(iBATCE_ASPC_comb, file = paste0(iBATCE_ASPC_prefix,"Seurat.rds"))
