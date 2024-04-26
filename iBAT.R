source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================
GetCellComposition(iBAT_merged)
GetCellComposition(iBAT_Ad_merged)
GetCellComposition(iBAT_ASPC_merged)

table(iBAT_merged$orig.ident)
PlotObjMetrices(Seurat_Obj = iBAT_merged, OutPrefix= iBAT_prefix, ColorUse = iBAT_SAMcolors)
#===============================================================================
# 1. iBAT integration
#===============================================================================
iBAT_prefix <- "iBAT_ACells_"
if(file.exists(paste0(iBAT_prefix,"Seurat.rds"))){
  iBAT_merged <- readRDS(paste0(iBAT_prefix,"Seurat.rds"))
} else {
  iBAT1mo_data <- ReadCB_h5("../_rawdata/iBAT-1mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iBAT2mo_data <- ReadCB_h5("../_rawdata/iBAT-2mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iBAT6mo_data <- ReadCB_h5("../_rawdata/iBAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iBAT1mo_Obj <- CreateSeuratObject(counts = iBAT1mo_data, project = "iBAT-1mo", min.cells = 10)
  iBAT2mo_Obj <- CreateSeuratObject(counts = iBAT2mo_data, project = "iBAT-2mo", min.cells = 10)
  iBAT6mo_Obj <- CreateSeuratObject(counts = iBAT6mo_data, project = "iBAT-6mo", min.cells = 10)
  iBAT_merged <- RunSeurat(ObjList = c(iBAT1mo_Obj,  iBAT2mo_Obj, iBAT6mo_Obj), OutPrefix = iBAT_prefix, NormalizationMethod = "VST", 
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15, DoubletPercent = 0.08)
  
  Idents(iBAT_merged) <- iBAT_merged$integrated_snn_res.0.1
  DimPlot(iBAT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = iBAT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_merged, ColorUse = hughie_color, OutPrefix = iBAT_prefix)
  saveRDS(iBAT_merged, file = paste0(iBAT_prefix,"Seurat.rds"))
}

iBAT_merged <- RenameIdents(iBAT_merged,  `0` = "Adipocyte", `1` = "Endothelial Cell", `2` = "Adipocyte",`3` = "ASPC", 
                              `4` = "Adipocyte", `5` = "Pericyte", `6` = "Smooth Muscle Cell",`7` = "Macrophage")
iBAT_merged$Myannotation <- iBAT_merged@active.ident
iBAT_cellRank <- c("ASPC","Adipocyte","Macrophage","Endothelial Cell","Pericyte","Smooth Muscle Cell")
iBAT_colors <- c("#0072B5FF","#BC3C29FF","#7E6148FF","#631879FF","#3B4992FF","#F39B7FFF")
names(iBAT_colors) <- iBAT_cellRank
iBAT_SAMcolors <- c( "#A2D48E","#47A93A", "#729E4E")
Idents(iBAT_merged) <- factor(Idents(iBAT_merged), levels = iBAT_cellRank)
iBAT_merged$orig.ident <- factor(iBAT_merged$orig.ident, levels = c("iBAT-1mo", "iBAT-2mo", "iBAT-6mo"))

PlotGenes <- c("Pdgfra","Nova1", # ASPC
               "Plin1","Adipoq","Acvr1c", # Adipocyte
               "Dock2","Mrc1",# Macrophage
               "Cyyr1","Flt1", # Endothelial Cell
               "Rgs5","Abcc9", # Pericyte
               "Mylk4","Trdn"# Smooth Muscle Cell
)

(p1 <- VlnPlot(iBAT_merged, features = PlotGenes, fill.by = "ident", cols = iBAT_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 10, y.text = 10) + NoLegend())
(p2 <- VlnPlot(iBAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iBAT_SAMcolors, assay = "RNA") + FontSize(x.text = 10, y.text = 10))
p1 + p2
ggsave(paste0(iBAT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 12, height = 7)

(p1 <- DimPlot(iBAT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = iBAT_colors, label.color = "white"))
(p2 <- DimPlot(iBAT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iBAT_SAMcolors))
p1 + p2
ggsave(paste0(iBAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(iBAT_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T,  cols = iBAT_colors, label.color = "white", order = T)
ggsave(paste0(iBAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iBAT_merged, OutPrefix = iBAT_prefix, ColorUse1 = iBAT_colors, ColorUse2 = iBAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_merged, ColorUse = iBAT_colors, OutPrefix = iBAT_prefix, NMarkers = 15)
DEG_enrichment(Seruat_DEG_file = paste0(iBAT_prefix, "Allmarkers_", Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4,cellRanks = iBAT_cellRank)

saveRDS(iBAT_merged, file = paste0(iBAT_prefix,"Seurat.rds"))

#===============================================================================
# 2. iBAT-Adipocytes
#===============================================================================
iBAT_Ad_prefix <- "iBAT_Adipocytes_"
if(file.exists(paste0(iBAT_Ad_prefix,"Seurat.rds"))){
  iBAT_Ad_merged <- readRDS(paste0(iBAT_Ad_prefix,"Seurat.rds"))
} else {
  iBAT_Ad_merged <- ScaleData(subset(iBAT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(iBAT_Ad_merged)
  ggsave(paste0(iBAT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(iBAT_Ad_merged) <- iBAT_Ad_merged$integrated_snn_res.0.15
  DimPlot(iBAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iBAT_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBAT_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_Ad_merged, ColorUse = hughie_color, OutPrefix = iBAT_Ad_prefix, NMarkers = 20)
}

iBAT_Ad_merged <- RenameIdents(iBAT_Ad_merged, `0` = "iBAT-ad1", `1` = "iBAT-ad3", `2` = "iBAT-ad2", `3` = "iBAT-ad4")
iBAT_Ad_merged$Myannotation_Ad <- iBAT_Ad_merged@active.ident
iBAT_AdRanks <- c("iBAT-ad1","iBAT-ad2","iBAT-ad3","iBAT-ad4")
iBAT_Ad_colors <- c("#BC3C29FF", "#EE4C97FF","#FB9A99","#33A02C")
names(iBAT_Ad_colors) <- iBAT_AdRanks
Idents(iBAT_Ad_merged) <- factor(Idents(iBAT_Ad_merged), levels = iBAT_AdRanks)

PlotGenes <- c("Plin1","Ucp1","Cidea","Cyp2e1","Slit3","Slc7a10","Aldh1a1","Lep")
(p1 <- VlnPlot(iBAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = iBAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iBAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iBAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iBAT_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(iBAT_Ad_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = iBAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(iBAT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iBAT_SAMcolors))
p1 + p2
ggsave(paste0(iBAT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(iBAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = iBAT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(iBAT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iBAT_Ad_merged, OutPrefix = iBAT_Ad_prefix, ColorUse1 = iBAT_Ad_colors, ColorUse2 = iBAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_Ad_merged, ColorUse = iBAT_Ad_colors, OutPrefix = iBAT_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iBAT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4)

FeaturePlot(iBAT_Ad_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(iBAT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 8, height = 7.5)
saveRDS(iBAT_Ad_merged, file = paste0(iBAT_Ad_prefix,"Seurat.rds"))

Calculate_moduleScore(Seurat_Obj = subset(iBAT_Ad_merged, subset = orig.ident == "iBAT-2mo"), ColorUse = iBAT_Ad_colors, NtopGene = 50, Ylimits = c(0, 1),
                      Seruat_DEG_file = "../1. rATM_126mo/rATM_Adipocytes_2mo_Allmarkers.csv", OutPrefix = paste0(iBAT_Ad_prefix, "2mo_mPRAT_"))

Calculate_moduleScore(Seurat_Obj = subset(iBAT_Ad_merged, subset = orig.ident == "iBAT-2mo"), ColorUse = iBAT_Ad_colors, NtopGene = 50, Ylimits = c(0, 1),
                      Seruat_DEG_file = "../2. Compare/2020Sun_nature_Allmarkers.csv", OutPrefix = paste0(iBAT_Ad_prefix, "2020Sun_"))

#------------------- Trajectory analysis ---------------------------
cds <- as.cell_data_set(iBAT_Ad_merged,assay = "RNA")
rowData(cds)$gene_short_name <- rownames(rowData(cds))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds), verbose = T)

pdf(paste0(iBAT_Ad_prefix,"Sudotime_",Sys.Date(),".pdf"), width = 5, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, cell_size = 0.2, label_roots = T, label_leaves=F, 
           label_branch_points=F, graph_label_size=1, trajectory_graph_segment_size = 0.5, trajectory_graph_color = "black")
dev.off()

#===============================================================================
# 3. iBAT-ASPC
#===============================================================================
iBAT_ASPC_prefix <- "iBAT_ASPC_"
if(file.exists(paste0(iBAT_ASPC_prefix,"Seurat.rds"))){
  iBAT_ASPC_merged <- readRDS(paste0(iBAT_ASPC_prefix,"Seurat.rds"))
} else {
  iBAT_ASPC_merged <- ScaleData(subset(iBAT_merged, idents = c("ASPC"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(iBAT_ASPC_merged)
  ggsave(paste0(iBAT_ASPC_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(iBAT_ASPC_merged) <- iBAT_ASPC_merged$integrated_snn_res.0.15
  DimPlot(iBAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = iBAT_ASPC_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = iBAT_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_ASPC_merged, ColorUse = hughie_color, OutPrefix = iBAT_ASPC_prefix)
}

iBAT_ASPC_merged <- RenameIdents(iBAT_ASPC_merged, `0` = "iBAT-aspc2", `1` = "iBAT-aspc3", `2` = "iBAT-aspc1")
iBAT_ASPC_merged$Myannotation_ASPC <- iBAT_ASPC_merged@active.ident
iBAT_ASPCRanks <- c("iBAT-aspc1","iBAT-aspc2","iBAT-aspc3")
iBAT_ASPC_colors <- c("#1F78B4","#FDBF6F","#CAB2D6")
names(iBAT_ASPC_colors) <- cellRank
Idents(iBAT_ASPC_merged) <- factor(Idents(iBAT_ASPC_merged), levels = iBAT_ASPCRanks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3", "Pdgfrb","Fmo2","Pparg","Cd36","Cidea","Ucp1")

(p1 <- VlnPlot(iBAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = iBAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(iBAT_ASPC_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = iBAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(iBAT_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(iBAT_ASPC_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = iBAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(iBAT_ASPC_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = iBAT_SAMcolors))
p1 + p2
ggsave(paste0(iBAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(iBAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = iBAT_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(iBAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = iBAT_ASPC_merged, OutPrefix = iBAT_ASPC_prefix, ColorUse1 = iBAT_ASPC_colors, ColorUse2 = iBAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = iBAT_ASPC_merged, ColorUse = iBAT_ASPC_colors, OutPrefix = iBAT_ASPC_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(iBAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = iBAT_ASPCRanks)

FeaturePlot(iBAT_ASPC_merged, features =  c("Pdgfra","Dpp4","Pdgfrb","Cidea"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(iBAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(iBAT_ASPC_merged, file = paste0(iBAT_ASPC_prefix,"Seurat.rds"))
