source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================

GetCellComposition(lPRAT_merged)
GetCellComposition(lPRAT_Ad_merged)
GetCellComposition(lPRAT_ASPC_merged)
GetCellComposition(lPRAT_Mac_merged)

#===============================================================================
# 1. lPRAT integration
#===============================================================================
lPRAT_prefix <- "lPRAT_ACells_"
if(file.exists(paste0(lPRAT_prefix,"Seurat.rds"))){
  lPRAT_merged <- readRDS(paste0(lPRAT_prefix,"Seurat.rds"))
} else {
  lPRAT2mo_data <- ReadCB_h5("../_rawdata/lPRAT-2mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  lPRAT6mo_data <- ReadCB_h5("../_rawdata/lPRAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  lPRAT2mo_Obj <- CreateSeuratObject(counts = lPRAT2mo_data, project = "lPRAT-2mo", min.cells = 10)
  lPRAT6mo_Obj <- CreateSeuratObject(counts = lPRAT6mo_data, project = "lPRAT-6mo", min.cells = 10)
  lPRAT_merged <- RunSeurat(ObjList = c(lPRAT2mo_Obj, lPRAT6mo_Obj), OutPrefix = lPRAT_prefix, NormalizationMethod = "VST", 
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)
  Idents(lPRAT_merged) <- lPRAT_merged$integrated_snn_res.0.15
  DimPlot(lPRAT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = lPRAT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = lPRAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_merged, ColorUse = hughie_color, OutPrefix = lPRAT_prefix)
}

lPRAT_merged <- RenameIdents(lPRAT_merged,`0` = "Adipocyte", `1` = "Mesothelial Cell",`2` = "Macrophage",`3` = "ASPC", 
                              `4` = "Adipocyte",`5` = "Macrophage",`6`="Mesothelial Cell",`7`="Macrophage")
lPRAT_merged$Myannotation <- lPRAT_merged@active.ident
lPRAT_cellRank <- c("ASPC","Adipocyte","Macrophage","Mesothelial Cell")
lPRAT_colors <- c("#0072B5FF","#BC3C29FF","#7E6148FF","#631879FF")
names(lPRAT_colors) <- lPRAT_cellRank
lPRAT_SAMcolors <- c("#A6CEE3", "#1F78B4")
Idents(lPRAT_merged) <- factor(Idents(lPRAT_merged), levels= lPRAT_cellRank)
lPRAT_merged$orig.ident <- factor(lPRAT_merged$orig.ident, levels = c("lPRAT-2mo","lPRAT-6mo"))

PlotGenes <- c("Pdgfra","Nova1", # ASPC
               "Plin1","Acvr1c","Lep","Cidea","Ucp1", # Adipocyte
               "Dock2","Mrc1",# Macrophage
               "Muc16","Msln"# Mesothelial Cell
)

(p1 <- VlnPlot(lPRAT_merged, features = PlotGenes, fill.by = "ident", cols = lPRAT_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(lPRAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = lPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(lPRAT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 10, height = 8)

(p1 <- DimPlot(lPRAT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T, cols = lPRAT_colors, label.color = "white"))
(p2 <- DimPlot(lPRAT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = lPRAT_SAMcolors))
p1 + p2
ggsave(paste0(lPRAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(lPRAT_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, cols = lPRAT_colors, label.color = "white", order = T)
ggsave(paste0(lPRAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 6)

Plot_Cell_compoistion(Seurat_Obj = lPRAT_merged, OutPrefix = lPRAT_prefix, ColorUse1 = lPRAT_colors, ColorUse2 = lPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_merged, ColorUse = lPRAT_colors, OutPrefix = lPRAT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(lPRAT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = lPRAT_cellRank)
PlotObjMetrices(Seurat_Obj = lPRAT_merged, OutPrefix= lPRAT_prefix, ColorUse = lPRAT_SAMcolors)

saveRDS(lPRAT_merged, file = paste0(lPRAT_prefix,"Seurat.rds"))

#===============================================================================
# 2. lPRAT-Adipocytes
#===============================================================================
lPRAT_Ad_prefix <- "lPRAT_Adipocytes_"
if(file.exists(paste0(lPRAT_Ad_prefix,"Seurat.rds"))){
  lPRAT_Ad_merged <- readRDS(paste0(lPRAT_Ad_prefix,"Seurat.rds"))
} else {
  lPRAT_Ad_merged <- ScaleData(subset(lPRAT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(lPRAT_Ad_merged)
  ggsave(paste0(lPRAT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  saveRDS(lPRAT_Ad_merged, file = paste0(lPRAT_Ad_prefix,"Seurat.rds"))
  
  Idents(lPRAT_Ad_merged) <- lPRAT_Ad_merged$integrated_snn_res.0.2
  DimPlot(lPRAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = lPRAT_Ad_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = lPRAT_Ad_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_Ad_merged, ColorUse = hughie_color, OutPrefix = lPRAT_Ad_prefix, NMarkers = 40)
}

lPRAT_Ad_merged <- RenameIdents(lPRAT_Ad_merged, `0` = "lPRAT-ad2", `1` = "lPRAT-ad1", `2` = "lPRAT-ad3", `3` = "lPRAT-ad4")
lPRAT_Ad_merged$Myannotation_Ad <- lPRAT_Ad_merged@active.ident
lPRAT_AdRanks <- c("lPRAT-ad1","lPRAT-ad2","lPRAT-ad3","lPRAT-ad4")
lPRAT_Ad_colors <- c( "#729E4E","#FDBF6F","#A20056FF","#EE4C97FF")
names(lPRAT_Ad_colors) <- lPRAT_AdRanks
Idents(lPRAT_Ad_merged) <- factor(Idents(lPRAT_Ad_merged), levels = lPRAT_AdRanks)

PlotGenes <- c("Adrb3","Lpl","Cd36", #Lipolytic
               "Acly","Fasn",#Lipogenesis
               "Anxa2","Ighm","Fgf1","Hif1a", #endocytosis and lipid transport
               "Irf2","Gulp1","Mitf",
               "Cmss1","Lars2","Camk1d","Cdk8"
)

(p1 <- VlnPlot(lPRAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = lPRAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(lPRAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = lPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(lPRAT_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(lPRAT_Ad_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = lPRAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(lPRAT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = lPRAT_SAMcolors))
p1 + p2
ggsave(paste0(lPRAT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(lPRAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = lPRAT_Ad_colors, label.color = "white", order = T)
ggsave(paste0(lPRAT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = lPRAT_Ad_merged, OutPrefix = lPRAT_Ad_prefix, ColorUse1 = lPRAT_Ad_colors, ColorUse2 = lPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_Ad_merged, ColorUse = lPRAT_Ad_colors, OutPrefix = lPRAT_Ad_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(lPRAT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4,cellRanks = lPRAT_AdRanks)

FeaturePlot(lPRAT_Ad_merged, features =  c("Lpl","Acly","Anxa2","Cdk8"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(lPRAT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)

##------------------- Merge with Mandurp adipocytes --------------------------- 
AdMerged_prefix <- "Merge2021Mandrup_Ad_"
if(file.exists(paste0(Merged_prefix,"Seurat.rds"))){
  AdMerged_merged <- readRDS(paste0(AdMerged_prefix,"Seurat.rds"))
} else {
  MandrupAd <- readRDS(paste0("2021Mandrup_Seurat.rds"))
  ObjList = c(lPRAT_Ad_merged, MandrupAd)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  AdMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=5))
  AdMerged_merged$orig.ident <- gsub("-[2|6]mo","",AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- gsub("LFD","eWAT",AdMerged_merged$orig.ident)
  AdMerged_merged$orig.ident <- as.factor(AdMerged_merged$orig.ident)
  saveRDS(AdMerged_merged, file = paste0(AdMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Plin1","Ucp1","Cidea","Adipoq","Cyp2e1","Aldh1a1",
               "Pparg","Acaca","Elovl6","Acly","Igf2r",
               "Nnat","Lrp3","Car3","Abcd2",
               "Hif1a","Nedd9","Gadd45g","Lep","Rab7")

Idents(AdMerged_merged) <- AdMerged_merged$integrated_snn_res.0.2
(p00 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p01 <- DimPlot(AdMerged_merged, reduction = "umap", group.by = "orig.ident", label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p02 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = c("#DC0000FF","#1F78B4"), stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))
(p00 + p01)/p02
ggsave(paste0(AdMerged_prefix,"UMAP-mergedIdentity_",Sys.Date(),".pdf"), width = 12, height = 12)
Plot_Cell_compoistion(Seurat_Obj = AdMerged_merged, OutPrefix = Merged_prefix, ColorUse1 = rev(hughie_color[1:4]), ColorUse2 = c("#DC0000FF","#1F78B4"))

Idents(AdMerged_merged) <- AdMerged_merged$Subtype
(p10 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#FF7F00","#6A3D9A","#20854EFF")))
(pb <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA")) 

Idents(AdMerged_merged) <- AdMerged_merged$Myannotation_Ad
Idents(AdMerged_merged) <- factor(Idents(AdMerged_merged), levels = lPRAT_AdRanks)
(p20 <- DimPlot(AdMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = lPRAT_Ad_colors))
(pc <- VlnPlot(AdMerged_merged, features = PlotGenes, fill.by = "ident",cols = lPRAT_Ad_colors, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p01)/(p10)/(p20)
ggsave(paste0(AdMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

#===============================================================================
# 3. lPRAT-ASPC
#===============================================================================
lPRAT_ASPC_prefix <- "lPRAT_ASPC_"
if(file.exists(paste0(lPRAT_ASPC_prefix,"Seurat.rds"))){
  lPRAT_ASPC_merged <- readRDS(paste0(lPRAT_ASPC_prefix,"Seurat.rds"))
} else {
  lPRAT_ASPC_merged <- ScaleData(subset(lPRAT_merged, idents = c("ASPC"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  
  Idents(lPRAT_ASPC_merged) <- lPRAT_ASPC_merged$integrated_snn_res.0.35
  DimPlot(lPRAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = lPRAT_ASPC_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = lPRAT_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_ASPC_merged, ColorUse = hughie_color, OutPrefix = lPRAT_ASPC_prefix)
}

lPRAT_ASPC_merged <- RenameIdents(lPRAT_ASPC_merged, `0` = "lPRAT-aspc1",`1` = "lPRAT-aspc1", `2` = "lPRAT-aspc2")
lPRAT_ASPC_merged$Myannotation_ASPC <- lPRAT_ASPC_merged@active.ident
lPRAT_ASPCRanks <- c("lPRAT-aspc1","lPRAT-aspc2")
lPRAT_ASPC_colors <- c("#0072B5FF", "#E18727FF")
names(lPRAT_ASPC_colors) <- lPRAT_ASPCRanks
Idents(lPRAT_ASPC_merged) <- factor(Idents(lPRAT_ASPC_merged), levels = lPRAT_ASPCRanks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36")

(p1 <- VlnPlot(lPRAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = lPRAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(lPRAT_ASPC_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = lPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(lPRAT_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(lPRAT_ASPC_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = lPRAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(lPRAT_ASPC_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = lPRAT_SAMcolors))
p1 + p2
ggsave(paste0(lPRAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(lPRAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = lPRAT_ASPC_colors, label.color = "white", order = T)
ggsave(paste0(lPRAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = lPRAT_ASPC_merged, OutPrefix = lPRAT_ASPC_prefix, ColorUse1 = lPRAT_ASPC_colors, ColorUse2 = lPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_ASPC_merged, ColorUse = lPRAT_ASPC_colors, OutPrefix = lPRAT_ASPC_prefix,NMarkers = 50)
DEG_enrichment(Seruat_DEG_file = paste0(lPRAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4, cellRanks = lPRAT_ASPCRanks)

FeaturePlot(lPRAT_ASPC_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(lPRAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(lPRAT_ASPC_merged, file = paste0(lPRAT_ASPC_prefix,"Seurat.rds"))

##------------------- Merge with Mandurp ASPC --------------------------- 
ASPCMerged_prefix <- "Merge2021Mandrup_ASPC_"
if(file.exists(paste0(ASPCMerged_prefix,"Seurat.rds"))){
  ASPCMerged_merged <- readRDS(paste0(ASPCMerged_prefix,"Seurat.rds"))
} else {
  MandrupASPC <- readRDS(paste0("../../Public reference datasets/2021 Mandrup eWAT snRNAseq/eWAT_FAP.Rds"))
  MandrupASPC <- subset(MandrupASPC, subset = Dataset == "LFD_R1")
  MandrupASPC$orig.ident <- "eWAT"
  MandrupASPC$Subtype <- factor(MandrupASPC$Subtype, levels =c("FAP1","FAP2","FAP3","FAP4"))
  ObjList = c(lPRAT_ASPC_merged, MandrupASPC)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  ASPCMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 30) %>%  RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>%   FindClusters(resolution = seq(from=0, by=0.05, length=5))

  ASPCMerged_merged$orig.ident <- gsub("-[2|6]mo","", ASPCMerged_merged$orig.ident)
  ASPCMerged_merged$orig.ident <- factor(ASPCMerged_merged$orig.ident)
  saveRDS(ASPCMerged_merged, file = paste0(ASPCMerged_prefix,"Seurat.rds"))
}

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36")
Idents(ASPCMerged_merged) <- ASPCMerged_merged$integrated_snn_res.0.2
(p00 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p0 <- DimPlot(ASPCMerged_merged, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p1 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#BC3C29FF","#E18727FF","#0072B5FF","#20854EFF")))
(pa <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Subtype
Idents(ASPCMerged_merged) <- factor(Idents(ASPCMerged_merged) , levels =c("FAP1","FAP2","FAP3","FAP4"))
(p2 <- DimPlot(ASPCMerged_merged, reduction = "umap",split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#BC3C29FF","#E18727FF","#0072B5FF","#20854EFF")))
(pb <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(ASPCMerged_merged) <- ASPCMerged_merged$Myannotation_ASPC
(p3 <- DimPlot(ASPCMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#0072B5FF", "#E18727FF")))
(pc <- VlnPlot(ASPCMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p0)/(p2)/(p3)
ggsave(paste0(ASPCMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)

#===============================================================================
# 4. lPRAT-Macrophage
#===============================================================================
lPRAT_Mac_prefix <- "lPRAT_Mac_"
if(file.exists(paste0(lPRAT_Mac_prefix,"Seurat.rds"))){
  lPRAT_Mac_merged <- readRDS(paste0(lPRAT_Mac_prefix,"Seurat.rds"))
} else {
  lPRAT_Mac_merged <- ScaleData(subset(lPRAT_merged, idents = c("Macrophage"))) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(lPRAT_Mac_merged)
  ggsave(paste0(lPRAT_Mac_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  Idents(lPRAT_Mac_merged) <- lPRAT_Mac_merged$integrated_snn_res.0.1
  DimPlot(lPRAT_Mac_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = lPRAT_Mac_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = lPRAT_Mac_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_Mac_merged, ColorUse = hughie_color, OutPrefix = lPRAT_Mac_prefix)
}

lPRAT_Mac_merged <- RenameIdents(lPRAT_Mac_merged, `0` = "lPRAT-mac1",`1` = "lPRAT-mac2", `2` = "lPRAT-mac3")
lPRAT_Mac_merged$Myannotation_Macrophage <- lPRAT_Mac_merged@active.ident
lPRAT_MacRanks <- c("lPRAT-mac1","lPRAT-mac2","lPRAT-mac3")
lPRAT_Mac_colors <- c("#B15928","#F39B7FFF","#808180FF")
names(lPRAT_Mac_colors) <- lPRAT_MacRanks
Idents(lPRAT_Mac_merged) <- factor(Idents(lPRAT_Mac_merged), levels = lPRAT_MacRanks)

PlotGenes <- c("Adgre1","Mctp1","Trem2","Lpl","Apbb2","Atp6v0d2","Rbpj","Fn1","Cd163","Plekhg5")

(p1 <- VlnPlot(lPRAT_Mac_merged, features = PlotGenes, fill.by = "ident", cols = lPRAT_Mac_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(lPRAT_Mac_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = lPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(lPRAT_Mac_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(lPRAT_Mac_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = lPRAT_Mac_colors, label.color = "white"))
(p2 <- DimPlot(lPRAT_Mac_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = lPRAT_SAMcolors))
p1 + p2
ggsave(paste0(lPRAT_Mac_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(lPRAT_Mac_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = lPRAT_Mac_colors, label.color = "white", order = T)
ggsave(paste0(lPRAT_Mac_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 12, height = 5)

Plot_Cell_compoistion(Seurat_Obj = lPRAT_Mac_merged, OutPrefix = lPRAT_Mac_prefix, ColorUse1 = lPRAT_Mac_colors, ColorUse2 = lPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = lPRAT_Mac_merged, ColorUse = lPRAT_Mac_colors, OutPrefix = lPRAT_Mac_prefix)
DEG_enrichment(Seruat_DEG_file = paste0(lPRAT_Mac_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 15, filterLevel = 4)

FeaturePlot(lPRAT_Mac_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(lPRAT_Mac_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(lPRAT_Mac_merged, file = paste0(lPRAT_Mac_prefix,"Seurat.rds"))

##------------------- Merge with Mandurp Macrophage --------------------------- 
MacMerged_prefix <- "Merge2021Mandrup_Mac_"
if(file.exists(paste0(MacMerged_prefix,"Seurat.rds"))){
  MacMerged_merged <- readRDS(paste0(MacMerged_prefix,"Seurat.rds"))
} else {
  MandrupMacrophage <- readRDS(paste0("../../Public reference datasets/2021 Mandrup eWAT snRNAseq/eWAT_Annotated.Rds"))
  MandrupMacrophage <- subset(MandrupMacrophage, subset = Dataset == "LFD_R1")
  MandrupMacrophage <- subset(MandrupMacrophage, subset = Subtype %in% c("CEM","LAM","NPVM","P-LAM","PVM","RM"))
  MandrupMacrophage$orig.ident <- "eWAT"
  MandrupMacrophage$orig.ident <- as.factor(MandrupMacrophage$orig.ident)
  ObjList = c(lPRAT_Mac_merged, MandrupMacrophage)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  MacMerged_merged <- ScaleData(combined) %>% RunPCA(npcs = 30) %>%  RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>%   FindClusters(resolution = seq(from=0, by=0.05, length=5))

  MacMerged_merged$orig.ident <- gsub("-[2|6]mo","", MacMerged_merged$orig.ident)
  MacMerged_merged$orig.ident <- factor(MacMerged_merged$orig.ident)
  saveRDS(MacMerged_merged, file = paste0(MacMerged_prefix,"Seurat.rds"))
}

Idents(MacMerged_merged) <- MacMerged_merged$integrated_snn_res.0.1
(p00 <- DimPlot(MacMerged_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = hughie_color, label.color = "white"))
(p0 <- DimPlot(MacMerged_merged, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box = T, order = T, cols = c("#DC0000FF","#1F78B4")))
(p1 <- DimPlot(MacMerged_merged, reduction = "umap", group.by = "ident", ncol = 1, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pa <- VlnPlot(MacMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(MacMerged_merged) <- MacMerged_merged$Subtype
(p2 <- DimPlot(MacMerged_merged, reduction = "umap",split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(pb <- VlnPlot(MacMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

Idents(MacMerged_merged) <- MacMerged_merged$Myannotation_Macrophage
(p3 <- DimPlot(MacMerged_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box = T, order = T, cols = c("#B15928","#F39B7FFF","#808180FF")))
(pc <- VlnPlot(MacMerged_merged, features = PlotGenes, fill.by = "ident",cols = hughie_color, stack = T, split.by = "orig.ident",flip=T, assay = "RNA"))

(p00+p0)/(p2)/(p3)
ggsave(paste0(MacMerged_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 15)
