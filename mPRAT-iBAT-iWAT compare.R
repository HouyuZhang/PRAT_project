source("0.ScFunctions.R")
#===============================================================================
# 0. Basic stats
#===============================================================================
mPRATiBATiWAT_cells<- c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4", "iBAT-ad1", "iBAT-ad2", "iBAT-ad3","iBAT-ad4", "iWAT-ad1", "iWAT-ad2", "iWAT-ad3","iWAT-ad4","iWAT-ad5")
mPRATiBATiWAT_color <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF","#BC3C29FF", "#EE4C97FF","#0072B5FF","#20854EFF","#BC3C29FF","#FF7F00","#0072B5FF","#6A3D9A","#20854EFF")
names(mPRATiBATiWAT_color) <- mPRATiBATiWAT_cells
mPRATiBATiWAT_SAMcolors <- c("black","#E18727FF","#EE4C97FF")
PlotGenes <- c("Acvr1c","Ucp1","Cidea","Cyp2e1","Slit3","Slc7a10","Aldh1a1","Lep","Eya4")

#===============================================================================
# 1. mPRAT iBAT iWAT RT
#===============================================================================
mPRATiBATiWAT_RT_Ad_prefix <- "mPRATiBATiWAT_RT_Adipocyte_"
if(file.exists(paste0(mPRATiBATiWAT_RT_Ad_prefix, "Seurat.rds"))){
  mPRATiBATiWAT_RT_Ad_comb <- readRDS(paste0(mPRATiBATiWAT_RT_Ad_prefix, "Seurat.rds"))
} else {
  mPRAT_RT_Obj <- readRDS(paste0("mPRAT_CE_Adipocyte_Seurat.rds"))
  mPRAT_RT_Obj <- subset(mPRAT_RT_Obj, subset = orig.ident == "mPRAT-6mo-RT")
  mPRAT_RT_Obj$mPRATRT_anno <- mPRAT_RT_Obj$CEAnno
  
  iBAT_RT_Obj <- readRDS(paste0("iBAT_CE_Adipocyte_Seurat.rds"))
  iBAT_RT_Obj <- subset(iBAT_RT_Obj, subset = orig.ident == "iBAT-6mo-RT")
  iBAT_RT_Obj$iBATRT_anno <- iBAT_RT_Obj$Myad
  
  iWAT_RT_Obj <- readRDS(paste0("iWAT_RTvsCE_Adipocytes_Seurat.rds"))
  iWAT_RT_Obj <- subset(iWAT_RT_Obj, subset = orig.ident == "iWAT-RT")
  iWAT_RT_Obj$iWATRT_anno <- iWAT_RT_Obj$Myannotation_Ad
  
  ObjList <- c(mPRAT_RT_Obj, iWAT_RT_Obj,iBAT_RT_Obj)
  ObjList <- lapply(X = ObjList, FUN = function(Obj) {
    Obj <- NormalizeData(Obj, normalization.method = "LogNormalize", scale.factor = 10000)
    Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000)
    Obj <- ScaleData(Obj)
  })
  
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  mPRATiBATiWAT_RT_Ad_comb <- ScaleData(combined) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length = 10))
  clustree(mPRATiBATiWAT_RT_Ad_comb)
  ggsave(paste0(mPRATiBATiWAT_RT_Ad_prefix,"clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  mPRATiBATiWAT_RT_Ad_comb$orig.ident <- factor(mPRATiBATiWAT_RT_Ad_comb$orig.ident, levels = c("mPRAT-6mo-RT","iBAT-6mo-RT","iWAT-RT"))
  mPRATiBATiWAT_RT_Ad_comb$AllID <- mPRATiBATiWAT_RT_Ad_comb$mPRATRT_anno
  mPRATiBATiWAT_RT_Ad_comb@meta.data[mPRATiBATiWAT_RT_Ad_comb$orig.ident == "iWAT-RT",]$AllID <- mPRATiBATiWAT_RT_Ad_comb$iWATRT_anno[mPRATiBATiWAT_RT_Ad_comb$orig.ident == "iWAT-RT"]
  mPRATiBATiWAT_RT_Ad_comb@meta.data[mPRATiBATiWAT_RT_Ad_comb$orig.ident == "iBAT-6mo-RT",]$AllID <- mPRATiBATiWAT_RT_Ad_comb$iBATRT_anno[mPRATiBATiWAT_RT_Ad_comb$orig.ident == "iBAT-6mo-RT"]
  mPRATiBATiWAT_RT_Ad_comb$AllID <- factor(mPRATiBATiWAT_RT_Ad_comb$AllID, levels = mPRATiBATiWAT_cells)
  Idents(mPRATiBATiWAT_RT_Ad_comb) <- mPRATiBATiWAT_RT_Ad_comb$AllID
  saveRDS(mPRATiBATiWAT_RT_Ad_comb, file = paste0(mPRATiBATiWAT_RT_Ad_prefix,"Seurat.rds"))
}

(p1 <- DimPlot(mPRATiBATiWAT_RT_Ad_comb, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = mPRATiBATiWAT_SAMcolors))
(p2 <- DimPlot(mPRATiBATiWAT_RT_Ad_comb, reduction = "umap", group.by = "integrated_snn_res.0.4", label = F, repel = T, label.box = T, order = T, cols = hughie_color))
(p3 <- DimPlot(mPRATiBATiWAT_RT_Ad_comb, reduction = "umap", split.by = "orig.ident",group.by = "integrated_snn_res.0.4", ncol =3, label = F, order = T, cols = hughie_color))
(p4 <- DimPlot(mPRATiBATiWAT_RT_Ad_comb, reduction = "umap", split.by = "orig.ident", group.by = "AllID",ncol = 3, label = F, repel = T, label.box =T, order = T, cols = mPRATiBATiWAT_color))
(p1|p2)/p3/p4
ggsave(paste0(mPRATiBATiWAT_RT_Ad_prefix,"Compare",Sys.Date(),".pdf"), width = 14, height = 18)

VlnPlot(mPRATiBATiWAT_RT_Ad_comb, features = PlotGenes, fill.by = "ident", cols = mPRATiBATiWAT_color, stack = T, flip=T, assay = "RNA")+ NoLegend()
ggsave(paste0(mPRATiBATiWAT_RT_Ad_prefix,"Compare_VlnPlot_", Sys.Date(), ".pdf"), width = 6, height = 8)

MarkersPlot(Seurat_Obj = mPRATiBATiWAT_RT_Ad_comb, ColorUse = mPRATiBATiWAT_color, OutPrefix = mPRATiBATiWAT_RT_Ad_prefix, NMarkers = 20)

#===============================================================================
# 2. mPRAT iBAT iWAT CE
#===============================================================================
mPRATiBATiWAT_CE_Ad_prefix <- "mPRATiBATiWAT_CE_Adipocyte_"
if(file.exists(paste0(mPRATiBATiWAT_CE_Ad_prefix, "Seurat.rds"))){
  mPRATiBATiWAT_CE_Ad_comb <- readRDS(paste0(mPRATiBATiWAT_CE_Ad_prefix, "Seurat.rds"))
} else {
  mPRAT_CE_Obj <- readRDS(paste0("mPRAT_CE_Adipocyte_Seurat.rds"))
  mPRAT_CE_Obj <- subset(mPRAT_CE_Obj, subset = orig.ident == "mPRAT-6mo-CE")
  mPRAT_CE_Obj$mPRATCE_anno <- mPRAT_CE_Obj$CEAnno
  
  iBAT_CE_Obj <- readRDS(paste0("iBAT_CE_Adipocyte_Seurat.rds"))
  iBAT_CE_Obj <- subset(iBAT_CE_Obj, subset = orig.ident == "iBAT-6mo-CE")
  iBAT_CE_Obj$iBATCE_anno <- iBAT_CE_Obj$Myad
  
  iWAT_CE_Obj <- readRDS(paste0("iWAT_RT+CE/iWAT_RTvsCE_Adipocytes_Seurat.rds"))
  iWAT_CE_Obj <- subset(iWAT_CE_Obj, subset = orig.ident == "iWAT-CE")
  iWAT_CE_Obj$iWATCE_anno <- iWAT_CE_Obj$Myannotation_Ad
  
  ObjList <- c(mPRAT_CE_Obj, iWAT_CE_Obj,iBAT_CE_Obj)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  mPRATiBATiWAT_CE_Ad_comb <- ScaleData(combined) %>% RunPCA(npcs =30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(mPRATiBATiWAT_CE_Ad_comb)
  ggsave(paste0(mPRATiBATiWAT_CE_Ad_prefix,"clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  
  mPRATiBATiWAT_CE_Ad_comb$orig.ident <- factor(mPRATiBATiWAT_CE_Ad_comb$orig.ident, levels = c("mPRAT-6mo-CE","iBAT-6mo-CE","iWAT-CE"))
  mPRATiBATiWAT_CE_Ad_comb$AllID <- mPRATiBATiWAT_CE_Ad_comb$mPRATCE_anno
  mPRATiBATiWAT_CE_Ad_comb@meta.data[mPRATiBATiWAT_CE_Ad_comb$orig.ident == "iWAT-CE",]$AllID <- mPRATiBATiWAT_CE_Ad_comb$iWATCE_anno [mPRATiBATiWAT_CE_Ad_comb$orig.ident == "iWAT-CE"]
  mPRATiBATiWAT_CE_Ad_comb@meta.data[mPRATiBATiWAT_CE_Ad_comb$orig.ident == "iBAT-6mo-CE",]$AllID <- mPRATiBATiWAT_CE_Ad_comb$iBATCE_anno [mPRATiBATiWAT_CE_Ad_comb$orig.ident == "iBAT-6mo-CE"]
  mPRATiBATiWAT_CE_Ad_comb$AllID <- factor(mPRATiBATiWAT_CE_Ad_comb$AllID, levels = mPRATiBATiWAT_cells)
  Idents(mPRATiBATiWAT_CE_Ad_comb) <- mPRATiBATiWAT_CE_Ad_comb$AllID
  
  saveRDS(mPRATiBATiWAT_CE_Ad_comb, file = paste0(mPRATiBATiWAT_CE_Ad_prefix,"Seurat.rds"))
}

(p1 <- DimPlot(mPRATiBATiWAT_CE_Ad_comb, reduction = "umap", group.by = "orig.ident", ncol = 1, label = F, repel = T, label.box =T, order = T, cols = mPRATiBATiWAT_SAMcolors))
(p2 <- DimPlot(mPRATiBATiWAT_CE_Ad_comb, reduction = "umap", group.by = "integrated_snn_res.0.4",label = F, repel = T, label.box =T, order = T, cols = hughie_color))
(p3 <- DimPlot(mPRATiBATiWAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", group.by = "integrated_snn_res.0.4",ncol =3, label = F,  order = T, cols = hughie_color))
(p4 <- DimPlot(mPRATiBATiWAT_CE_Ad_comb, reduction = "umap", split.by = "orig.ident", group.by = "AllID",ncol = 3, label = F, repel = T, label.box =T, order = T, cols = mPRATiBATiWAT_color))
(p1|p2)/p3/p4
ggsave(paste0(mPRATiBATiWAT_CE_Ad_prefix,"Compare",Sys.Date(),".pdf"), width = 14, height = 18)

VlnPlot(mPRATiBATiWAT_CE_Ad_comb, features = PlotGenes, fill.by = "ident", cols = mPRATiBATiWAT_color, stack = T, flip=T, assay = "RNA")+ NoLegend()
ggsave(paste0(mPRATiBATiWAT_CE_Ad_prefix,"Compare_VlnPlot_", Sys.Date(), ".pdf"), width = 6, height = 8)
MarkersPlot(Seurat_Obj = mPRATiBATiWAT_CE_Ad_comb, ColorUse = mPRATiBATiWAT_color, OutPrefix = mPRATiBATiWAT_CE_Ad_prefix, NMarkers = 20)

#===============================================================================
# 3. compare overall identity
#===============================================================================
mPRATiBATiWAT_Ad_prefix <- "mPRATiBATiWAT_Adipocyte_"

Idents(mPRATiBATiWAT_CE_Ad_comb) <- mPRATiBATiWAT_CE_Ad_comb$MergedAnno
expression_averge_CE <- AverageExpression(mPRATiBATiWAT_CE_Ad_comb)
markers_CE <- read.csv("mPRATiBATiWAT_CE_Adipocyte_Allmarkers.csv") %>% filter(p_val_adj < 0.05)

CEmtx <- expression_averge_CE$RNA[unique(markers_CE$gene),-c(3)]
CEmtx <- t(scale(t(CEmtx))) %>% as.data.frame()
(pCE <- Heatmap(CEmtx, cluster_columns = T, cluster_rows = T,
                clustering_method_columns = "complete", clustering_method_rows = "complete", 
                clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean",
                show_row_names = F, show_column_names = T, column_names_rot = 90, column_names_centered = F, 
                column_names_side = "top", row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 11, fontfCE = "bold"),
                column_dend_height = unit(3, "cm"), column_dend_gp = gpar(lwd = 2),show_row_dend = F,
                col = colorRamp2(breaks = c(-2,0,2),colors = c("blue","white","red")),raster_device = "png",border = T,
                column_split = 2,heatmap_legend_param = list(title = "Scaled \ngene \nexpression", border = "black", legend_height = unit(4, "cm"))
))
pdf(paste0(mPRATiBATiWAT_Ad_prefix,"DEG_heatmap_CE_", Sys.Date(),".pdf"), width = 6, height = 10)
draw(pCE)
dev.off()

#===============================================================================
# 4. Cell-cell communication analysis
#===============================================================================
mPRATiBATiWAT_cells<- c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4", "iBAT-ad1", "iBAT-ad2", "iBAT-ad3","iBAT-ad4", "iWAT-ad1", "iWAT-ad2", "iWAT-ad3","iWAT-ad4","iWAT-ad5")
mPRATiBATiWAT_color <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF","#BC3C29FF", "#EE4C97FF","#0072B5FF","#20854EFF","#BC3C29FF","#FF7F00","#0072B5FF","#6A3D9A","#20854EFF")
names(mPRATiBATiWAT_color) <- mPRATiBATiWAT_cells

RT_meta <- mPRATiBATiWAT_RT_Ad_comb@meta.data

mPRAT_RT_meta <- subset(RT_meta, orig.ident == "mPRAT-6mo-RT")
mPRAT_RT_meta$Anno <- factor(mPRAT_RT_meta$mPRATRT_anno, levels = c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4"))
mPRAT_RT_data <- as.matrix(mPRATiBATiWAT_RT_Ad_comb@assays$RNA@data[,as.character(row.names(mPRAT_RT_meta))]) 
Run_CellChat(SeuratObj = mPRAT_RT_data, Seurat_MetaInfo = mPRAT_RT_meta, OutPrefix = "CellChat_mPRAT_RT", ColorUse = mPRATiBATiWAT_color[1:4])

iBAT_RT_meta <-  subset(RT_meta, orig.ident == "iBAT-6mo-RT")
iBAT_RT_meta$Anno <- factor(iBAT_RT_meta$iBATRT_anno, levels = c("iBAT-ad1", "iBAT-ad2", "iBAT-ad3","iBAT-ad4"))
iBAT_RT_data <- as.matrix(mPRATiBATiWAT_RT_Ad_comb@assays$RNA@data[,as.character(row.names(iBAT_RT_meta))]) 
Run_CellChat(SeuratObj = iBAT_RT_data, Seurat_MetaInfo = iBAT_RT_meta, OutPrefix = "CellChat_iBAT_RT", ColorUse = mPRATiBATiWAT_color[5:8])

iWAT_RT_meta <-  subset(RT_meta, orig.ident == "iWAT-RT")
iWAT_RT_meta$Anno <- factor(iWAT_RT_meta$iWATRT_anno, levels = c("iWAT-ad1", "iWAT-ad2", "iWAT-ad3","iWAT-ad4","iWAT-ad5"))
# iWAT_RT_meta <- iWAT_RT_meta[!iWAT_RT_meta$Anno %in% c("iWAT-ad1", "iWAT-ad2"),]
# iWAT_RT_meta$Anno <- droplevels(iWAT_RT_meta$Anno)
iWAT_RT_data <- as.matrix(mPRATiBATiWAT_RT_Ad_comb@assays$RNA@data[,as.character(row.names(iWAT_RT_meta))]) 
Run_CellChat(SeuratObj = iWAT_RT_data, Seurat_MetaInfo = iWAT_RT_meta, OutPrefix = "CellChat_iWAT_RT2", ColorUse = mPRATiBATiWAT_color[9:13])

CE_meta <- mPRATiBATiWAT_CE_Ad_comb@meta.data

mPRAT_CE_meta <- subset(CE_meta, orig.ident == "mPRAT-6mo-CE")
mPRAT_CE_meta$Anno <- factor(mPRAT_CE_meta$mPRATCE_anno, levels = c("mPRAT-ad1","mPRAT-ad2", "mPRAT-ad3","mPRAT-ad4"))
# mPRAT_CE_meta <- mPRAT_CE_meta[!mPRAT_CE_meta$Anno %in% c( "mPRAT-ad2"),]
# mPRAT_CE_meta$Anno <- droplevels(mPRAT_CE_meta$Anno)
mPRAT_CE_data <- as.matrix(mPRATiBATiWAT_CE_Ad_comb@assays$RNA@data[,as.character(row.names(mPRAT_CE_meta))]) 
Run_CellChat(SeuratObj = mPRAT_CE_data, Seurat_MetaInfo = mPRAT_CE_meta, OutPrefix = "CellChat_mPRAT_CE2", ColorUse = mPRATiBATiWAT_color[c(1,2,3,4)])

iBAT_CE_meta <-  subset(CE_meta, orig.ident == "iBAT-6mo-CE")
iBAT_CE_meta$Anno <- factor(iBAT_CE_meta$iBATCE_anno, levels = c("iBAT-ad1", "iBAT-ad2", "iBAT-ad3","iBAT-ad4"))
iBAT_CE_data <- as.matrix(mPRATiBATiWAT_CE_Ad_comb@assays$RNA@data[,as.character(row.names(iBAT_CE_meta))]) 
Run_CellChat(SeuratObj = iBAT_CE_data, Seurat_MetaInfo = iBAT_CE_meta, OutPrefix = "CellChat_iBAT_CE", ColorUse = mPRATiBATiWAT_color[5:8])

iWAT_CE_meta <-  subset(CE_meta, orig.ident == "iWAT-CE")
iWAT_CE_meta$Anno <- factor(iWAT_CE_meta$iWATCE_anno, levels = c("iWAT-ad1", "iWAT-ad2", "iWAT-ad3","iWAT-ad4","iWAT-ad5"))
iWAT_CE_data <- as.matrix(mPRATiBATiWAT_CE_Ad_comb@assays$RNA@data[,as.character(row.names(iWAT_CE_meta))]) 
Run_CellChat(SeuratObj = iWAT_CE_data, Seurat_MetaInfo = iWAT_CE_meta, OutPrefix = "CellChat_iWAT_CE", ColorUse = mPRATiBATiWAT_color[9:13])

mPRAT_list <- list(mPRAT_RT = readRDS("CellChat_mPRAT_RT.rds"), mPRAT_CE = readRDS("CellChat_mPRAT_CE2.rds"))
iBAT_list <- list(iBAT_RT = readRDS("CellChat_iBAT_RT.rds"), iBAT_CE = readRDS("CellChat_iBAT_CE.rds"))
iWAT_list <- list(iWAT_RT = readRDS("CellChat_iWAT_RT2.rds"), iWAT_CE = readRDS("CellChat_iWAT_CE.rds"))
All_list <- c(unlist(mPRAT_list), unlist(iBAT_list), unlist(iWAT_list))

Run_CellChat_Compare(CellChat_Obj_list = mPRAT_list, OutPrefix = "mPRAT_Compare", ColorUse = mPRATiBATiWAT_color[1:4])
Run_CellChat_Compare(CellChat_Obj_list = iBAT_list, OutPrefix = "iBAT_Compare", ColorUse = mPRATiBATiWAT_color[5:8])
Run_CellChat_Compare(CellChat_Obj_list = iWAT_list, OutPrefix = "iWAT_Compare", ColorUse = mPRATiBATiWAT_color[9:13])
Run_CellChat_Compare(CellChat_Obj_list = All_list, OutPrefix = "All_Compare", ColorUse = mPRATiBATiWAT_color)

CellChat_Obj_list = mPRAT_list; OutPrefix = "mPRAT_Compare"; ColorUse = mPRATiBATiWAT_color[1:4]
CellChat_Obj_list = iBAT_list; OutPrefix = "iBAT_Compare"; ColorUse = mPRATiBATiWAT_color[5:8]
CellChat_Obj_list = iWAT_list; OutPrefix = "iWAT_Compare"; ColorUse = mPRATiBATiWAT_color[9:13]
CellChat_Obj_list = All_list; OutPrefix = "All_Compare"; ColorUse = mPRATiBATiWAT_color

SigPathways <- c()
for (i in 1:length(CellChat_Obj_list)){SigPathways <- unique(c(SigPathways,CellChat_Obj_list[[i]]@netP$pathways))}
cellchat <- mergeCellChat(CellChat_Obj_list, add.names = names(CellChat_Obj_list),cell.prefix = TRUE)
pdf(paste0(OutPrefix, "_netVisual_circle_",Sys.Date(),".pdf"), width = 14, height = 4)
weight.max <- getMaxWeight(CellChat_Obj_list, attribute = c("idents","count"))
par(mfrow = c(1,6), xpd=TRUE)
for (i in 1:length(CellChat_Obj_list)) {
  if(i<3){ColorUse = mPRATiBATiWAT_color[1:4]
  } else if(i<5 & i>2) {ColorUse = mPRATiBATiWAT_color[5:8]
  } else{ColorUse = mPRATiBATiWAT_color[9:13]} 
  netVisual_circle(CellChat_Obj_list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2],vertex.label.cex = 1,edge.label.cex = 1.2,
                   edge.width.max = 5, title.name = paste0("Number of interactions - ", names(CellChat_Obj_list)[i]), color.use = ColorUse)
}
dev.off()
