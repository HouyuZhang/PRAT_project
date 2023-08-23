source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================
PlotObjMetrices(Seurat_Obj = subset(mPRAT_CE_merged, subset = orig.ident == "mPRAT-6mo-ACE"), OutPrefix = mPRAT_CE_prefix, ColorUse = mPRAT_CE_SMAcolors)
GetCellComposition(mPRAT_CE_merged)

#===============================================================================
# 1. mPRAT RT+ACE
#===============================================================================
mPRAT_CE_prefix <- "mPRAT_CE_ACells_"
if(file.exists(paste0(mPRAT_CE_prefix,"Seurat.rds"))){
  mPRAT_CE_merged <- readRDS(paste0(mPRAT_CE_prefix,"Seurat.rds"))
} else {
  mPRAT6mo_data <-    ReadCB_h5("../_rawdata/mPRAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = T)
  mPRAT6mo_CE_data <- ReadCB_h5("../_rawdata/mPRAT-6mo-ACE/CellBender_feature_bc_matrix_filtered.h5", use.names = T)
  mPRAT6mo_Obj <- CreateSeuratObject(counts = mPRAT6mo_data, project = "mPRAT-6mo-RT", min.cells = 10)
  mPRAT6mo_CE_Obj <- CreateSeuratObject(counts = mPRAT6mo_CE_data, project = "mPRAT-6mo-ACE", min.cells = 10)
  mPRAT_CE_merged <- RunSeurat(ObjList = c(mPRAT6mo_Obj, mPRAT6mo_CE_Obj), OutPrefix = mPRAT_CE_prefix,
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)
 
  Idents(mPRAT_CE_merged) <- mPRAT_CE_merged$integrated_snn_res.0.1
  DimPlot(mPRAT_CE_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = T, repel = T, label.box =F, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRAT_CE_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRAT_CE_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_merged, ColorUse = hughie_color, OutPrefix = mPRAT_CE_prefix)
}

mPRAT_CE_merged <- RenameIdents(mPRAT_CE_merged,  `0` = "Adipocyte", `1` = "Adipocyte", `2` = "ASPC",`3` = "Macrophage", `4` = "Endothelial Cell", `5` = "Mesothelial Cell", `6` = "Pericyte")
mPRAT_CE_celltypeRanks <- c("ASPC","Adipocyte","Macrophage","Endothelial Cell","Pericyte","Mesothelial Cell")
mPRAT_CE_colorMaps <- c("#0072B5FF","#BC3C29FF","#7E6148FF","#E18727FF","#EE4C97FF", "#631879FF")
names(mPRAT_CE_colorMaps) <- mPRAT_CE_celltypeRanks
mPRAT_CE_merged$Myannotation <- mPRAT_CE_merged@active.ident
Idents(mPRAT_CE_merged) <- factor(Idents(mPRAT_CE_merged), levels = mPRAT_CE_celltypeRanks)
mPRAT_CE_merged$orig.ident <- factor(mPRAT_CE_merged$orig.ident, levels = c("mPRAT-6mo-RT", "mPRAT-6mo-ACE"))

PlotGenes <- c("Pdgfra","Nova1", # ASPC
               "Plin1","Adipoq","Acvr1c",# Adipocyte
               "Ptprc","Mrc1",# Macrophage
               "Cyyr1","Flt1", # Endothelial Cell
               "Rgs5","Abcc9", # Pericyte
               "Muc16","Msln" # Mesothelial Cell
)
(p1 <- VlnPlot(mPRAT_CE_merged, features = PlotGenes, fill.by = "ident", cols = mPRAT_CE_colorMaps, stack = T, flip=T,assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRAT_CE_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = hughie_color, assay = "RNA")) + FontSize(x.text = 16, y.text = 16)
p1 + p2
ggsave(paste0(mPRAT_CE_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 20, height = 12)

(p1 <- DimPlot(mPRAT_CE_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mPRAT_CE_colorMaps, label.color = "white"))
(p2 <- DimPlot(mPRAT_CE_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = hughie_color))
p1 + p2
ggsave(paste0(mPRAT_CE_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 15, height = 6)

DimPlot(mPRAT_CE_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRAT_CE_colorMaps, label.color = "white", order = T)
ggsave(paste0(mPRAT_CE_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 13, height = 6)

Plot_Cell_compoistion(mPRAT_CE_merged, OutPrefix = mPRAT_CE_prefix, ColorUse1 = mPRAT_CE_colorMaps, ColorUse2 = mPRAT_CE_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_merged, ColorUse = mPRAT_CE_colorMaps, OutPrefix = mPRAT_CE_prefix, NMarkers = 20)

DEG_enrichment(Seruat_DEG_file = paste0(mPRAT_CE_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = mPRAT_CE_celltypeRanks)
saveRDS(mPRAT_CE_merged, file = paste0(mPRAT_CE_prefix,"Seurat.rds"))

#===============================================================================
# 2. mPRAT RT+ACE Adipocytes
#===============================================================================
mPRAT_CE_Ad_prefix <- "mPRAT_CE_Adipocyte_"
if(file.exists(paste0(mPRAT_CE_Ad_prefix, "Seurat.rds"))){
  mPRAT_CE_Ad_comb <- readRDS(paste0(mPRAT_CE_Ad_prefix, "Seurat.rds"))
} else {
  mPRAT_CE_Ad_comb <- subset(mPRAT_CE_merged, idents = c("Adipocyte")) %>% ScaleData() %>% RunPCA(npcs = 30) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>% FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=20))
  clustree(mPRAT_CE_Ad_comb)
  ggsave(paste0(mPRAT_CE_Ad_prefix, "clustree_", Sys.Date(), ".pdf"), width = 12, height = 10)
  
  Idents(mPRAT_CE_Ad_comb) <- mPRAT_CE_Ad_comb$integrated_snn_res.0.35
  DimPlot(mPRAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRAT_CE_Ad_comb, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRAT_CE_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_Ad_comb, ColorUse = hughie_color, OutPrefix = mPRAT_CE_Ad_prefix)
}

mPRAT_CE_Ad_comb <- RenameIdents(mPRAT_CE_Ad_comb, `0` = "mPRAT-ad1", `1` = "mPRAT-ad3", `2` = "mPRAT-ad3",`3` = "mPRAT-ad1",
                                 `4` = "mPRAT-ad4", `5` = "mPRAT-ad1", `6` = "mPRAT-ad2",`7` = "mPRAT-ad1",`8` = "mPRAT-ad3")
mPRAT_CE_Ad_celltypeRanks <- c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4")
mPRAT_CE_Ad_comb$ACEAnno <- mPRAT_CE_Ad_comb@active.ident
mPRAT_CE_Ad_colors <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")
names(mPRAT_CE_Ad_colors) <- mPRAT_CE_Ad_celltypeRanks
Idents(mPRAT_CE_Ad_comb) <- factor(Idents(mPRAT_CE_Ad_comb), levels = mPRAT_CE_Ad_celltypeRanks)
mPRAT_CE_Ad_comb$orig.ident <- factor(mPRAT_CE_Ad_comb$orig.ident, levels = c("mPRAT-6mo-RT", "mPRAT-6mo-ACE"))

PlotGenes <- c("Plin1","Ucp1","Cidea","Pank1","Gk","Cyp2e1","Tshr","Slc7a10","Aldh1a1")
(p1 <- VlnPlot(mPRAT_CE_Ad_comb, features = PlotGenes, fill.by = "ident", cols = mPRAT_CE_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRAT_CE_Ad_comb, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = hughie_color, assay = "RNA",pt.size = 1) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRAT_CE_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 16, height = 12)

(p1 <- DimPlot(mPRAT_CE_Ad_comb, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T,  cols = mPRAT_CE_Ad_colors, label.color = "white"))
(p2 <- DimPlot(mPRAT_CE_Ad_comb, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = hughie_color))
p1 + p2
ggsave(paste0(mPRAT_CE_Ad_prefix, "UMAP-ByClusters_", Sys.Date(), ".pdf"), width = 16, height = 8)

DimPlot(mPRAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRAT_CE_Ad_colors, label.color = "white", order = T)
ggsave(paste0(mPRAT_CE_Ad_prefix,"UMAP-Byorig.ident_", Sys.Date(), ".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = mPRAT_CE_Ad_comb, OutPrefix = mPRAT_CE_Ad_prefix, ColorUse1 = mPRAT_CE_Ad_colors, ColorUse2 = mPRAT_CE_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_Ad_comb, ColorUse = mPRAT_CE_Ad_colors, OutPrefix = mPRAT_CE_Ad_prefix, NMarkers = 30)

FeaturePlot(mPRAT_CE_Ad_comb, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"),ncol = 2, split.by = "orig.ident",min.cutoff = "q10", order = T, pt.size=0.3) 
ggsave(paste0(mPRAT_CE_Ad_prefix,"FeaturePlot_", Sys.Date(), ".pdf"), width = 10, height = 18)
saveRDS(mPRAT_CE_Ad_comb, file = paste0(mPRAT_CE_Ad_prefix,"Seurat.rds"))

#------------------- Compare annotation with time-series dataset ---------------------------
mPRAT_Ad_merged1 <- readRDS(paste0("../1. mPRAT_126mo/mPRAT_Adipocytes_Seurat.rds"))
TS_metadata <- subset(mPRAT_Ad_merged1, subset = orig.ident == "mPRAT-6mo")@meta.data
rownames(TS_metadata) <- gsub("_3","",rownames(TS_metadata))
mPRAT_Ad_combRT <- subset(mPRAT_CE_Ad_comb, subset = orig.ident == "mPRAT-6mo-RT")
Cells <- gsub("_1","",rownames(mPRAT_Ad_combRT@meta.data))
mPRAT_Ad_combRT@meta.data$TsAnno <- TS_metadata[Cells,]$Myannotation_Ad

(p1 <- DimPlot(mPRAT_Ad_combRT, reduction = "umap", group.by = "TsAnno", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = mPRAT_CE_Ad_colors))
(p2 <- DimPlot(mPRAT_Ad_combRT, reduction = "umap",label = F, repel = T, label.box =T, order = T, cols = mPRAT_CE_Ad_colors))
p1+p2
ggsave(paste0(mPRAT_CE_Ad_prefix,"Compare_Anno",Sys.Date(),".pdf"), width = 14, height = 6)

mPRAT_Ad_combRT@meta.data[,c("TsAnno","FinalID")] %>% dplyr::count(TsAnno, FinalID) %>% filter(TsAnno != "NA") %>% 
  ggplot(aes(axis1 = TsAnno, axis2 = FinalID, y = n)) +
  geom_alluvium(aes(fill = TsAnno),alpha=0.5) + geom_stratum(aes(fill = TsAnno)) + geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")) + theme_void() + theme(legend.position = "none")
ggsave(paste0(mPRAT_CE_Ad_prefix,"SankeyPlot_",Sys.Date(),".pdf"), width = 4, height = 8)

#------------------- DEG ---------------------------
mPRAT_CE_Ad_comb$CT_SAMPLE <- paste(mPRAT_CE_Ad_comb$orig.ident, Idents(mPRAT_CE_Ad_comb), sep = "_")
Idents(mPRAT_CE_Ad_comb) <- "CT_SAMPLE"
mito.genes <- grep(pattern = "^mt-|^Hb|^EN", rownames(mPRAT_CE_Ad_comb@assays$RNA), value = TRUE)
mPRAT_CE_Ad_comb <- subset(mPRAT_CE_Ad_comb, features = setdiff(rownames(mPRAT_CE_Ad_comb@assays$RNA), mito.genes))

#List to store all pair-wise DEGs
DEGmat <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(DEGmat) <-c("gene","p_val", "avg_log2FC", "pct.1", "pct.2","p_val_adj","change","cellType")

for (ad in c("ad1","ad3","ad4")){
  Seurat_obj_DEG <- FindMarkers(mPRAT_CE_Ad_comb, ident.1 = paste0("mPRAT-6mo-RT_mPRAT-",ad), ident.2 = paste0("mPRAT-6mo-ACE_mPRAT-", ad), logfc.threshold = 0, min.pct = 0.1, assay = "RNA")
  Seurat_obj_DEG$change <- "Stable"
  Seurat_obj_DEG$change[Seurat_obj_DEG$p_val_adj >= 0.05] <- "Nosig"
  Seurat_obj_DEG$change[Seurat_obj_DEG$p_val_adj < 0.05 & Seurat_obj_DEG$avg_log2FC >= 0.25] <- "SigDown"
  Seurat_obj_DEG$change[Seurat_obj_DEG$p_val_adj < 0.05 & Seurat_obj_DEG$avg_log2FC < -0.25] <- "Sigup"
  Seurat_obj_DEG$cellType <- ad
  Seurat_obj_DEG <- Seurat_obj_DEG %>% rownames_to_column("gene")
  DEGmat <- rbind(DEGmat, Seurat_obj_DEG)
}
table(DEGmat$cellType, DEGmat$change) %>% t()
write.csv(DEGmat, file = paste0(mPRAT_CE_Ad_prefix,"DEGs_list_",Sys.Date(),".csv"), quote = F)
DEGmat <- read.csv("mPRAT_CE_Adipocyte2_DEGs_list.csv", row.names = 1)

top.marker.tmp <- DEGmat %>% filter(p_val_adj < 0.05) %>% group_by(cellType)
top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = 10, order_by = avg_log2FC)
top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = 10, order_by = avg_log2FC)
top.marker <- rbind(top.marker.max,top.marker.min)

DEGmat %>% filter(change != "Nosig") %>% 
ggplot(aes(x = cellType, y=avg_log2FC)) +
  geom_jitter(aes(color = change), size=0.5) + 
  ggrepel::geom_text_repel(data = top.marker, ggplot2::aes(x = cellType, y = avg_log2FC, label = gene, size = 1), max.overlaps = 50) +
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF", "grey")) +
  theme_classic() + 
  labs(x="", y="Average Log2 fold change (FC)") +
  theme(
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.text.x=element_blank(),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.text.y = element_text(color="black", size=16, face="bold")
  )
ggsave(paste0(mPRAT_CE_Ad_prefix,"DEG_scatterPlot_",Sys.Date(),".pdf"), width = 10, height = 6)

DEG_list <- list()
DEG_list[["ad1_up"]] <- DEGmat %>% filter(change == "Sigup" & cellType == "ad1") %>% pull(gene)
DEG_list[["ad1_down"]] <- DEGmat %>% filter(change == "SigDown" & cellType == "ad1") %>% pull(gene)
DEG_list[["ad3_up"]] <- DEGmat %>% filter(change == "Sigup" & cellType == "ad3") %>% pull(gene)
DEG_list[["ad3_down"]] <- DEGmat %>% filter(change == "SigDown" & cellType == "ad3") %>% pull(gene)
DEG_list[["ad4_up"]] <- DEGmat %>% filter(change == "Sigup" & cellType == "ad4") %>% pull(gene)
DEG_list[["ad4_down"]] <- DEGmat %>% filter(change == "SigDown" & cellType == "ad4") %>% pull(gene)
write.csv(do.call(cbind, DEG_list), file = "DEGs_list.csv", quote = F)

#===============================================================================
# 3. mPRAT RT ACE ASPC
#===============================================================================
mPRAT_CE_ASPC_prefix <- "mPRAT_CE_ASPC_"
if(file.exists(paste0(mPRAT_CE_ASPC_prefix,"Seurat.rds"))){
  mPRAT_CE_ASPC_combined <- readRDS(paste0(mPRAT_CE_ASPC_prefix,"Seurat.rds"))
} else {
  mPRAT_CE_ASPC_combined <- ScaleData(subset(mPRAT_CE_merged, idents = c("ASPC"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  clustree(mPRAT_CE_ASPC_combined)
  ggsave(paste0(mPRAT_CE_ASPC_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)

  Idents(mPRAT_CE_ASPC_combined) <- mPRAT_CE_ASPC_combined$integrated_snn_res.0.15
  DimPlot(mPRAT_CE_ASPC_combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRAT_CE_ASPC_combined, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRAT_CE_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_ASPC_combined, ColorUse = hughie_color, OutPrefix = mPRAT_CE_ASPC_prefix)
}

mPRAT_CE_ASPC_combined <- RenameIdents(mPRAT_CE_ASPC_combined, `0` = "mPRAT-aspc2", `1` = "mPRAT-aspc1", `2` = "mPRAT-aspc1", `3` = "mPRAT-aspc3")
mPRAT_CE_ASPC_combined$Myannotation_ASPC <- mPRAT_CE_ASPC_combined@active.ident
mPRAT_CE_ASPCRanks <- c("mPRAT-aspc1","mPRAT-aspc2","mPRAT-aspc3")
mPRAT_CE_ASPC_colors <- c("#0072B5FF","#E18727FF", "#6A3D9A")
names(mPRAT_CE_ASPC_colors) <- mPRAT_CE_ASPCRanks
Idents(mPRAT_CE_ASPC_combined) <- factor(Idents(mPRAT_CE_ASPC_combined), levels = mPRAT_CE_ASPCRanks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36","Cidea")

(p1 <- VlnPlot(mPRAT_CE_ASPC_combined, features = PlotGenes, fill.by = "ident", cols = mPRAT_CE_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRAT_CE_ASPC_combined, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRAT_CE_SAMcolors, assay = "RNA", split.plot = T) + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRAT_CE_ASPC_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 16, height = 12)

(p1 <- DimPlot(mPRAT_CE_ASPC_combined, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mPRAT_CE_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(mPRAT_CE_ASPC_combined, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRAT_CE_SAMcolors))
p1 + p2
ggsave(paste0(mPRAT_CE_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 16, height = 8)

DimPlot(mPRAT_CE_ASPC_combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRAT_CE_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(mPRAT_CE_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = mPRAT_CE_ASPC_combined, OutPrefix = mPRAT_CE_ASPC_prefix, ColorUse1 = mPRAT_CE_ASPC_colors, ColorUse2 = mPRAT_CE_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_CE_ASPC_combined, ColorUse = mPRAT_CE_ASPC_colors, OutPrefix = mPRAT_CE_ASPC_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRAT_CE_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4,cellRanks = mPRAT_CE_ASPCRanks)

FeaturePlot(mPRAT_CE_ASPC_combined, features =  c("Pdgfra","Dpp4","Pdgfrb","Fabp4"), min.cutoff = "q9", order = T, pt.size=0.3)
ggsave(paste0(mPRAT_CE_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(mPRAT_CE_ASPC_combined, file = paste0(mPRAT_CE_ASPC_prefix,"Seurat.rds"))
