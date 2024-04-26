source("0.ScFunctions.R")
#===============================================================================
# 1. mPRAT-iWAT RT Adipocytes
#===============================================================================
mPRATiWATRT_prefix <- "mPRAT-iWAT_RT_"
if(file.exists(paste0(mPRATiWATRT_prefix,"Seurat.rds"))){
  mPRATiWATRT_merged <- readRDS(paste0(mPRATiWATRT_prefix,"Seurat.rds"))
} else {
  mPRAT_data <- ReadCB_h5("../_rawdata/mPRAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  iWAT_data <- ReadCB_h5("../_rawdata/iWAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  mPRAT_Obj <- CreateSeuratObject(counts = mPRAT_data, project = "mPRAT-RT", min.cells = 10)
  iWAT_Obj <- CreateSeuratObject(counts = iWAT_data, project = "iWAT-RT", min.cells = 10)
  
  mPRATiWATRT_merged <- RunSeurat(ObjList = c(mPRAT_Obj, iWAT_Obj), Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), 
                                 OutPrefix = mPRATiWATRT_prefix, NormalizationMethod = "VST", MitoPercent = 15)

  Idents(mPRATiWATRT_merged) <- mPRATiWATRT_merged$integrated_snn_res.0.05
  DimPlot(mPRATiWATRT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = mPRATiWATRT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRATiWATRT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWATRT_merged, ColorUse = hughie_color, OutPrefix = mPRATiWATRT_prefix)
}

mPRATiWATRT_merged <- RenameIdents(mPRATiWATRT_merged,  `0` = "Adipocyte", `1` = "B Cell", `2` = "T Cell",`3` = "Macrophage", `4` = "ASPC", 
                                    `5` = "Endothelial Cell", `6` = "B Cell",`7` = "T Cell", `8` = "Pericyte", `9` = "Kidney-related Cell",
                                    `10` ="Neutrophill",`11` ="Dendritic Cell")
mPRATiWATRT_merged$Myannotation <- mPRATiWATRT_merged@active.ident
mPRATiWATRT_cellRanks <- c("ASPC","Adipocyte","Macrophage","T Cell","B Cell","Neutrophill","Dendritic Cell","Endothelial Cell","Pericyte", "Kidney-related Cell")
mPRATiWATRT_color <- c("#0072B5FF","#BC3C29FF","#7E6148FF","#E18727FF","#EE4C97FF","#F39B7FFF","#B09C85FF","#631879FF","#3B4992FF", "#808180FF")
names(mPRATiWATRT_color) <- mPRATiWATRT_cellRanks
Idents(mPRATiWATRT_merged) <- factor(Idents(mPRATiWATRT_merged), levels= mPRATiWATRT_cellRanks)
mPRATiWATRT_merged$orig.ident <- factor(mPRATiWATRT_merged$orig.ident, levels = c("mPRAT-RT","iWAT-RT"))
mPRATiWATRT_SAMcolors <- c("#6A3D9A","#0072B5FF")

PlotGenes <- c("Pdgfra","Nova1", #ASPC
               "Plin1","Acvr1c",#Adipocytes
               "Dock2","Ptprc","Adgre1","Mrc1", #Macrophage
               "Prkcq","Themis",# T Cell
               "Pax5","Cr2", # B Cell
               "Csf3r", #Neutrophill
               "Flt3", #Dendritic Cell
               "Cyyr1","Flt1", #Endothelial
               "Abcc9","Rgs5" #Pericyte
)

(p1 <- VlnPlot(mPRATiWATRT_merged, features = PlotGenes, fill.by = "ident", cols = mPRATiWATRT_color, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRATiWATRT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRATiWATRT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRATiWATRT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width= 14, height = 8)

(p1 <- DimPlot(mPRATiWATRT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mPRATiWATRT_color, label.color = "white"))
(p2 <- DimPlot(mPRATiWATRT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRATiWATRT_SAMcolors))
p1 + p2
ggsave(paste0(mPRATiWATRT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mPRATiWATRT_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRATiWATRT_color, label.color = "white", order = T)
ggsave(paste0(mPRATiWATRT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 14, height = 12)

Plot_Cell_compoistion(Seurat_Obj = subset(mPRATiWATRT_merged, idents = mPRATiWATRT_cellRanks[c(1:9)]), OutPrefix = mPRATiWATRT_prefix, ColorUse1 = mPRATiWATRT_color, ColorUse2 = mPRATiWATRT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWATRT_merged, ColorUse = mPRATiWATRT_color, OutPrefix = mPRATiWATRT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRATiWATRT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = mPRATiWATRT_cellRanks)
saveRDS(mPRATiWATRT_merged, file = paste0(mPRATiWATRT_prefix,"Seurat.rds"))

#===============================================================================
# 2. mPRAT-iWAT RT Adipocytes
#===============================================================================
mPRATiWATRT_Ad_prefix <- "mPRAT-iWAT_RT_Adipocytes_"
if(file.exists(paste0(mPRATiWATRT_Ad_prefix,"Seurat.rds"))){
  mPRATiWATRT_Ad_merged <- readRDS(paste0(mPRATiWATRT_Ad_prefix,"Seurat.rds"))
} else {
  mPRATiWATRT_Ad_merged <- ScaleData(subset(mPRATiWATRT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  clustree(mPRATiWATRT_Ad_merged)
  ggsave(paste0(mPRATiWATRT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$integrated_snn_res.0.2
  DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRATiWATRT_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRATiWATRT_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWATRT_Ad_merged, ColorUse = hughie_color, OutPrefix = mPRATiWATRT_Ad_prefix)
}

mPRATiWATRT_Ad_merged <- RenameIdents(mPRATiWATRT_Ad_merged, `0` = "ad1", `1` = "ad3", `2` = "ad4",`3` = "ad2",`4` = "ad5")
mPRATiWATRT_merged$Myannotation_Ad <- mPRATiWATRT_Ad_merged@active.ident
mPRATiWATRT_AdRanks <- c("ad1","ad2","ad3","ad4","ad5")
mPRATiWATRT_Ad_colors <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF","#0072B5FF")
names(mPRATiWATRT_Ad_colors) <- mPRATiWATRT_AdRanks
Idents(mPRATiWATRT_Ad_merged) <- factor(Idents(mPRATiWATRT_Ad_merged), levels = mPRATiWATRT_AdRanks)

PlotGenes <- c("Acvr1c","Adrb3","Prdm16","Dio2","Ucp1","Cidea","Aldh1a1","Lep","Trhde","Atp2a2")

(p1 <- VlnPlot(mPRATiWATRT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mPRATiWATRT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRATiWATRT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRATiWATRT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRATiWATRT_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 16, height = 12)

# Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$Myannotation_Ad
(p1 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = mPRATiWATRT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRATiWATRT_SAMcolors))
p1 + p2
ggsave(paste0(mPRATiWATRT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 12, height = 5)

DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = mPRATiWATRT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(mPRATiWATRT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 8, height = 8)

Plot_Cell_compoistion(Seurat_Obj = mPRATiWATRT_Ad_merged, OutPrefix = mPRATiWATRT_Ad_prefix, ColorUse1 = mPRATiWATRT_Ad_colors, ColorUse2 = mPRATiWATRT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRATiWATRT_Ad_merged, ColorUse = Ad_colors, OutPrefix = mPRATiWATRT_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRATiWATRT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = mPRATiWATRT_AdRanks)

FeaturePlot(mPRATiWATRT_Ad_merged, features =  c("Ucp1","Cidea","Esrrg","Ar","Cyp2e1","Trhde"), min.cutoff = "q10", order = T, pt.size = 0.3)
ggsave(paste0(mPRATiWATRT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 8, height = 10)
saveRDS(mPRATiWATRT_Ad_merged, file = paste0(mPRATiWATRT_Ad_prefix,"Seurat.rds"))

#------------------- annotation by individual dataset ---------------------------
mPRAT_Ad <- readRDS(paste0("../1. mPRAT_126mo/mPRAT_Adipocytes_Seurat.rds"))
mPRAT_Ad_metadata <- subset(mPRAT_Ad, subset = orig.ident == "mPRAT-6mo")@meta.data
rownames(mPRAT_Ad_metadata) <- gsub("_3","_1",rownames(mPRAT_Ad_metadata))

iWAT_Ad <- readRDS(paste0("../7. iWAT_RT+CE/iWAT_RTvsACE_Adipocytes_Seurat.rds"))
iWAT_Ad_metadata <- subset(iWAT_Ad, subset = orig.ident == "iWAT-RT")@meta.data
rownames(iWAT_Ad_metadata) <- gsub("_1","_2",rownames(iWAT_Ad_metadata))

mPRATiWATRT_Ad_merged$mPRAT_iWAT_anno <- "NA"
over1 <- intersect(rownames(mPRATiWATRT_Ad_merged@meta.data),rownames(mPRAT_Ad_metadata))
mPRATiWATRT_Ad_merged$mPRAT_old_anno <- "NA"
mPRATiWATRT_Ad_merged@meta.data[over1,]$mPRAT_old_anno <- as.vector(mPRAT_Ad_metadata[over1,]$Myannotation_Ad)
mPRATiWATRT_Ad_merged@meta.data[over1,]$mPRAT_iWAT_anno <- as.vector(mPRAT_Ad_metadata[over1,]$Myannotation_Ad)
over2 <- intersect(rownames(mPRATiWATRT_Ad_merged@meta.data),rownames(iWAT_Ad_metadata))
mPRATiWATRT_Ad_merged$iWAT_old_anno <- "NA"
mPRATiWATRT_Ad_merged@meta.data[over2,]$iWAT_old_anno <- as.vector(iWAT_Ad_metadata[over2,]$Myannotation_Ad)
mPRATiWATRT_Ad_merged@meta.data[over2,]$mPRAT_iWAT_anno <- as.vector(iWAT_Ad_metadata[over2,]$Myannotation_Ad)

mPRATiWATRT_Ad_merged$orig.ident <- as.factor(mPRATiWATRT_Ad_merged$orig.ident)

Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$Myannotation
(p0 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = c("#EE4C97FF","#1B1919FF")))
(p1 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = hughie_color))

mPRATiWATRT_Ad_merged$mPRAT_old_anno <- factor(mPRATiWATRT_Ad_merged$mPRAT_old_anno, levels = c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4", "NA"))
Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$mPRAT_old_anno
(p2 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, 
               order = T, cols = c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF","grey")))

mPRATiWATRT_Ad_merged$iWAT_old_anno <- factor(mPRATiWATRT_Ad_merged$iWAT_old_anno, levels = c("iWAT-ad1", "iWAT-ad2", "iWAT-ad3", "iWAT-ad4", "iWAT-ad5", "NA"))
Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$iWAT_old_anno
(p3 <- DimPlot(mPRATiWATRT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, 
               order = T, cols = c("#BC3C29FF","#FF7F00","#0072B5FF","#6A3D9A","#20854EFF","grey")))
(p0|p1)/p2/p3
ggsave(paste0("Compare_RTACE",Sys.Date(),".pdf"), width = 11, height = 14)

mPRATiWATRT_Ad_merged$mPRAT_iWAT_anno <- factor(mPRATiWATRT_Ad_merged$mPRAT_iWAT_anno, 
                                                levels = c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4", "iWAT-ad1", "iWAT-ad2", "iWAT-ad3", "iWAT-ad4", "iWAT-ad5", "NA"))
Idents(mPRATiWATRT_Ad_merged) <- mPRATiWATRT_Ad_merged$mPRAT_iWAT_anno

PlotGenes <- c("Acvr1c","Ucp1","Cidea","Cyp2e1","Slit3","Slc7a10","Aldh1a1","Lep")
VlnPlot(mPRATiWATRT_Ad_merged, features = PlotGenes, fill.by = "ident", idents = c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4",  "iWAT-ad3", "iWAT-ad4", "iWAT-ad5"),
        cols = c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF","#0072B5FF","#6A3D9A","#20854EFF","grey"), stack = T, flip=T, assay = "RNA")+ NoLegend()
ggsave(paste0("Compare_RTACE_VlnPlot_", Sys.Date(), ".pdf"), width = 6, height = 8)

