source("0.ScFunctions.R")

#===============================================================================
# 0. Basic stats
#===============================================================================
GetCellComposition(mPRAT_merged)
GetCellComposition(mPRAT_Ad_merged)
GetCellComposition(mPRAT_ASPC_merged)

table(mPRAT_merged$orig.ident)
table(mPRAT_Ad_merged$orig.ident)
table(mPRAT_ASPC_merged$orig.ident)

PlotObjMetrices(Seurat_Obj = mPRAT_merged, OutPrefix= mPRAT_prefix, ColorUse = mPRAT_SAMcolors)

#===============================================================================
# 1. mPRAT integration
#===============================================================================
mPRAT_prefix <- "mPRAT_ACells_"
if(file.exists(paste0(mPRAT_prefix,"Seurat.rds"))){
  mPRAT_merged <- readRDS(paste0(mPRAT_prefix,"Seurat.rds"))
} else {
  mPRAT1mo_data <- ReadCB_h5("../_rawdata/mPRAT-1mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  mPRAT2mo_data <- ReadCB_h5("../_rawdata/mPRAT-2mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  mPRAT6mo_data <- ReadCB_h5("../_rawdata/mPRAT-6mo/CellBender_feature_bc_matrix_filtered.h5", use.names = TRUE)
  mPRAT1mo_Obj <- CreateSeuratObject(counts = mPRAT1mo_data, project = "mPRAT-1mo", min.cells = 10)
  mPRAT2mo_Obj <- CreateSeuratObject(counts = mPRAT2mo_data, project = "mPRAT-2mo", min.cells = 10)
  mPRAT6mo_Obj <- CreateSeuratObject(counts = mPRAT6mo_data, project = "mPRAT-6mo", min.cells = 10)
  mPRAT_merged <- RunSeurat(ObjList = c(mPRAT1mo_Obj, mPRAT2mo_Obj, mPRAT6mo_Obj), OutPrefix = mPRAT_prefix, NormalizationMethod = "VST", 
                             Filter_nFeature_RNA = c(500, 6000), Filter_nCount_RNA = c(1000, 50000), MitoPercent = 15)
  Idents(mPRAT_merged) <- mPRAT_merged$integrated_snn_res.0.1
  DimPlot(mPRAT_merged, reduction = "umap", ncol = 2, label = F, split.by = "orig.ident", repel = T, cols = hughie_color, label.box =T, order = T)
  Plot_Cell_compoistion(Seurat_Obj = mPRAT_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRAT_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_merged, ColorUse = hughie_color, OutPrefix = mPRAT_prefix)
}

mPRAT_merged <- RenameIdents(mPRAT_merged,  `0` = "Adipocyte", `1` = "Adipocyte", `2` = "ASPC",`3` = "Macrophage", `4` = "Adipocyte",
                             `5` = "Endothelial Cell",`6` = "Pericyte", `7` = "Mesothelial Cell")
mPRAT_merged$Myannotation <- mPRAT_merged@active.ident
mPRAT_cellRanks <- c("ASPC","Adipocyte","Macrophage","Endothelial Cell","Pericyte","Mesothelial Cell")
mPRAT_colorMaps <- c("#0072B5FF","#BC3C29FF","#7E6148FF","#E18727FF","#EE4C97FF", "#631879FF")
names(mPRAT_colorMaps) <- mPRAT_cellRanks
mPRAT_SAMcolors <- c("#C4ABD2","#9774B6", "#6A3D9A")
Idents(mPRAT_merged) <- factor(Idents(mPRAT_merged), levels = mPRAT_cellRanks)
mPRAT_merged$orig.ident <- factor(mPRAT_merged$orig.ident, levels = c("mPRAT-1mo", "mPRAT-2mo","mPRAT-6mo"))

PlotGenes <- c("Pdgfra","Nova1", # ASPC
               "Plin1","Adipoq","Acvr1c",# Adipocyte
               "Ptprc","Mrc1",# Macrophage
               "Cyyr1","Flt1", # Endothelial Cell
               "Rgs5","Abcc9", # Pericyte
               "Muc16","Msln" # Mesothelial Cell
)

(p1 <- VlnPlot(mPRAT_merged, features = PlotGenes, fill.by = "ident", cols = mPRAT_colorMaps, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 10, y.text = 10) + NoLegend())
(p2 <- VlnPlot(mPRAT_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 10, y.text = 10))
p1 + p2
ggsave(paste0(mPRAT_prefix, "Selected_markers_VlnPlot_",Sys.Date(),".pdf"), width = 12, height = 7)

(p1 <- DimPlot(mPRAT_merged, reduction = "umap", group.by = c("ident"), label = T, repel = T, label.box =T,  cols = mPRAT_colorMaps, label.color = "white"))
(p2 <- DimPlot(mPRAT_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRAT_SAMcolors))
p1 + p2
ggsave(paste0(mPRAT_prefix, "UMAP-ByClusters_",Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mPRAT_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, repel = T, label.box =T, cols = mPRAT_colorMaps, label.color = "white", order = T)
ggsave(paste0(mPRAT_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = subset(mPRAT_merged, idents = mPRAT_cellRanks), OutPrefix = mPRAT_prefix, ColorUse1 = mPRAT_colorMaps, ColorUse2 = mPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_merged, ColorUse = mPRAT_colorMaps, OutPrefix = mPRAT_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRAT_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = mPRAT_cellRanks)

#===============================================================================
# 2. mPRAT-Adipocytes
#===============================================================================
mPRAT_Ad_prefix <- "mPRAT_Adipocytes_"
if(file.exists(paste0(mPRAT_Ad_prefix,"Seurat.rds"))){
  mPRAT_Ad_merged <- readRDS(paste0(mPRAT_Ad_prefix,"Seurat.rds"))
} else {
  mPRAT_Ad_merged <- ScaleData(subset(mPRAT_merged, idents = c("Adipocyte"))) %>% RunPCA(npcs = 30) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(mPRAT_Ad_merged)
  ggsave(paste0(mPRAT_Ad_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
  Idents(mPRAT_Ad_merged) <- mPRAT_Ad_merged$integrated_snn_res.0.45
  DimPlot(mPRAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
}

mPRAT_Ad_merged <- RenameIdents(mPRAT_Ad_merged, `0` = "mPRAT-ad1", `1` = "mPRAT-ad2", `2` = "mPRAT-ad3", `3` = "mPRAT-ad2",
                                 `4` = "mPRAT-ad4",`5` = "mPRAT-ad3",`6` = "mPRAT-ad1",`7` = "mPRAT-ad1")
mPRAT_Ad_merged$Myannotation_Ad <- mPRAT_Ad_merged@active.ident
mPRAT_Ad_Ranks <- c("mPRAT-ad1","mPRAT-ad2","mPRAT-ad3","mPRAT-ad4")
mPRAT_Ad_colors <- c("#BC3C29FF","#FB9A99","#B2DF8A","#20854EFF")
names(mPRAT_Ad_colors) <- mPRAT_Ad_Ranks
Idents(mPRAT_Ad_merged) <- factor(Idents(mPRAT_Ad_merged), levels = mPRAT_Ad_Ranks)

PlotGenes <- c("Plin1","Acvr1c","Ucp1","Cidea","Pank1","Gk","Cpt1b",
               "Cyp2e1","Slit3","Tshr","Slc7a10","Aldh1a1","Lep","Npr3")

(p1 <- VlnPlot(mPRAT_Ad_merged, features = PlotGenes, fill.by = "ident", cols = mPRAT_Ad_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRAT_Ad_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRAT_Ad_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(mPRAT_Ad_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = mPRAT_Ad_colors, label.color = "white"))
(p2 <- DimPlot(mPRAT_Ad_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRAT_SAMcolors))
p1 + p2
ggsave(paste0(mPRAT_Ad_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mPRAT_Ad_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, label.box =T, cols = mPRAT_Ad_colors, label.color = "white", repel = T,order = T)
ggsave(paste0(mPRAT_Ad_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = mPRAT_Ad_merged, OutPrefix = mPRAT_Ad_prefix, ColorUse1 = mPRAT_Ad_colors, ColorUse2 = mPRAT_SAMcolors)
MarkersPlot(Seurat_Obj = mPRAT_Ad_merged, ColorUse = mPRAT_Ad_colors, OutPrefix = mPRAT_Ad_prefix, NMarkers = 40)

DEG_enrichment(Seruat_DEG_file = paste0(mPRAT_Ad_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = mPRAT_Ad_Ranks)
FeaturePlot(mPRAT_Ad_merged, features =  c("Ucp1","Cidea","Cyp2e1","Aldh1a1"), min.cutoff = "q10", order = F,pt.size=0.1, ncol = 2)
ggsave(paste0(mPRAT_Ad_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 8, height = 8)
saveRDS(mPRAT_Ad_merged, file = paste0(mPRAT_Ad_prefix,"Seurat.rds"))

#------------------- Trajectory analysis ---------------------------
cds <- as.cell_data_set(mPRAT_Ad_merged,assay ="RNA")
rowData(cds)$gene_short_name <- rownames(rowData(cds))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

pdf(paste0(mPRAT_Ad_prefix,"Sudotime_",Sys.Date(),".pdf"), width = 5, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=F, cell_size = 0.35, label_roots = F, label_leaves=F, 
           label_branch_points=F, graph_label_size=3, trajectory_graph_segment_size = 1, trajectory_graph_color = "black")
dev.off()

mPRAT_Ad_merged <- AddMetaData(object = mPRAT_Ad_merged, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "pseudotime")
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 20)
write.csv(modulated_genes, file = paste0(mPRAT_Ad_prefix, "modulated_genes_",Sys.Date(),".csv"), quote = F)

genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.1))
genes <- setdiff(genes, grep("^mt-",genes,value = T))
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes, rownames(rowData(cds))), order(pseudotime(cds))]
metainfor <- colData(cds)[order(pseudotime(cds)),]
metainfor$pseudotime <- pseudotime(cds)[order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

pdf(paste0(mPRAT_Ad_prefix,"Sudotime_heatmap_",Sys.Date(),".pdf"), width = 8, height = 8)
ha = columnAnnotation(CellType = metainfor$Myannotation_Ad,col = list(CellType = mPRAT_Ad_colors))
p <- Heatmap(pt.matrix, name  = "z-score", col = colorRamp2(seq(from=-2,to=2, length=11), rev(brewer.pal(11, "Spectral"))),
  show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 6), raster_device = "png",
  row_title_rot = 0, km = 3, cluster_rows = T, cluster_row_slices = FALSE, cluster_columns = FALSE, bottom_annotation = ha)
draw(p)
dev.off()
saveRDS(p, paste0(mPRAT_Ad_prefix,"Sudotime_heatmap_res_",Sys.Date(),".rds"))

#------------------- GO analysis ---------------------------
GeneNames <- list()
GeneNames[["decrease"]] <- rownames(pt.matrix)[row_order(p)[[2]]]
GeneNames[["midlle"]] <- rownames(pt.matrix)[row_order(p)[[3]]]
GeneNames[["increase"]] <- rownames(pt.matrix)[row_order(p)[[1]]]

compGO <- compareCluster(geneCluster = GeneNames, fun= "enrichGO", pvalueCutoff = 0.05, keyType = "SYMBOL", pAdjustMethod = "BH", OrgDb = org.Mm.eg.db, readable = T, ont = "ALL")
filterLevel = 4
compGO_filtered <- gofilter(compGO, level = filterLevel)
write.csv(compGO_filtered@compareClusterResult, file = paste0(mPRAT_Ad_prefix,"_pseudotime_GO-BP_list_",Sys.Date(),".csv"))

(p <- dotplot(compGO, showCategory = showCategoryNum, includeAll=T, split = "ONTOLOGY") + theme(axis.text.x = element_text(angle = 45, hjust = 1))) 
ggsave(paste0(mPRAT_Ad_prefix,"pseudotime_GO-BP_dotplot_slected_Filter",filterLevel,"_",Sys.Date(),".pdf"), width = 6, height = 7)

Go_BP_up <- enrichGO(gene = GeneNames[[1]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = 'BH', qvalueCutoff = 0.05)
Go_BP_down <- enrichGO(gene = GeneNames[[2]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = 'BH',  qvalueCutoff = 0.05)
GO_BP_All <- enrichGO(gene = GeneNames[[3]], OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = 'BH',  qvalueCutoff = 0.05)

p1 <- barplot(gofilter(Go_BP_up,level = filterLevel), showCategory = showCategoryNum,  title = paste0("GO-BP analysis of up DEG")) +  
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- barplot(gofilter(Go_BP_down,level = filterLevel), showCategory = showCategoryNum, title  = paste0("GO-BP analysis of down DEG"))+  
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- barplot(gofilter(GO_BP_All, level = filterLevel), showCategory = showCategoryNum, title  = paste0("GO-BP analysis"))+  
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p1/p2/p3)
ggsave(paste0(OutPrefix,"GO-BP analysis seperate",Sys.Date(),".pdf"),width = 14, height = 30)

#===============================================================================
# 3. mPRAT-ASPC
#===============================================================================
mPRAT_ASPC_prefix <- "mPRAT_ASPC_"
if(file.exists(paste0(mPRAT_ASPC_prefix,"Seurat.rds"))){
  mPRAT_ASPC_merged <- readRDS(paste0(mPRAT_ASPC_prefix,"Seurat.rds"))
} else {
  mPRAT_ASPC_merged <- ScaleData(subset(mPRAT_merged, idents = c("ASPC"))) %>% RunPCA(npcs = 30) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>% FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = seq(from=0, by=0.05, length=10))
  clustree(mPRAT_ASPC_merged)
  ggsave(paste0(mPRAT_ASPC_prefix, "clustree_",Sys.Date(),".pdf"), width = 12, height = 10)
 
  Idents(mPRAT_ASPC_merged) <- mPRAT_ASPC_merged$integrated_snn_res.0.1
  DimPlot(mPRAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 2, label = F, repel = T, label.box =T, order = T, cols = hughie_color)
  Plot_Cell_compoistion(Seurat_Obj = mPRAT_ASPC_merged, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = mPRAT_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_ASPC_merged, ColorUse = hughie_color, OutPrefix = mPRAT_ASPC_prefix)
}

mPRAT_ASPC_merged <- RenameIdents(mPRAT_ASPC_merged, `0` = "mPRAT-aspc1", `1` = "mPRAT-aspc2", `2` = "mPRAT-aspc3")
mPRAT_ASPC_merged$Myannotation_ASPC <- mPRAT_ASPC_merged@active.ident
mPRAT_ASPC_Ranks <- c("mPRAT-aspc1","mPRAT-aspc2","mPRAT-aspc3")
mPRAT_ASPC_colors <- c("#0072B5FF","#E18727FF", "#6A3D9A")
names(mPRAT_ASPC_colors) <- mPRAT_ASPC_Ranks
Idents(mPRAT_ASPC_merged) <- factor(Idents(mPRAT_ASPC_merged), levels = mPRAT_ASPC_Ranks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36","Cidea","Ucp1")

(p1 <- VlnPlot(mPRAT_ASPC_merged, features = PlotGenes, fill.by = "ident", cols = mPRAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(mPRAT_ASPC_merged, features = PlotGenes, stack=T, flip=T, split.by = "orig.ident", cols = mPRAT_SAMcolors, assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(mPRAT_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(mPRAT_ASPC_merged, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = mPRAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(mPRAT_ASPC_merged, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = mPRAT_SAMcolors))
p1 + p2
ggsave(paste0(mPRAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 14, height = 6)

DimPlot(mPRAT_ASPC_merged, reduction = "umap", split.by = "orig.ident", ncol = 3, label = F, label.box =T, cols = mPRAT_ASPC_colors, label.color = "white", repel = T,order = T)
ggsave(paste0(mPRAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 16, height = 6)

Plot_Cell_compoistion(Seurat_Obj = mPRAT_ASPC_merged, OutPrefix = mPRAT_ASPC_prefix, ColorUse1 = mPRAT_ASPC_colors, ColorUse2 = mPRAT_SAMcolors)
TopMarkersGenes <- MarkersPlot(Seurat_Obj = mPRAT_ASPC_merged, ColorUse = mPRAT_ASPC_colors, OutPrefix = mPRAT_ASPC_prefix, NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(mPRAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = mPRAT_ASPC_Ranks)

FeaturePlot(mPRAT_ASPC_merged, features = c("Pdgfra","Dpp4","Pdgfrb","Cidea"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(mPRAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(mPRAT_ASPC_merged, file = paste0(mPRAT_ASPC_prefix,"Seurat.rds"))

#===============================================================================
# 4. Merge mRPAT-ASPC with iBAT-ASPC
#===============================================================================
PRiBAT_ASPC_prefix <- "PRiBAT_ASPC_"
if(file.exists(paste0(PRiBAT_ASPC_prefix,"Seurat.rds"))){
  PRiBAT_ASPC_cob <- readRDS(paste0(PRiBAT_ASPC_prefix,"Seurat.rds"))
} else {
  ObjList = c(mPRAT_ASPC_merged, iBAT_ASPC_merged)
  features <- SelectIntegrationFeatures(object.list = ObjList, nfeatures = 2000)
  anchors <-  FindIntegrationAnchors(object.list = ObjList, anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  PRiBAT_ASPC_cob <- ScaleData(combined) %>% RunPCA(npcs = 30) %>% RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = seq(from=0, by=0.05, length=10))
  PRiBAT_ASPC_cob$orig.ident <- gsub("-[0-9]mo","",PRiBAT_ASPC_cob$orig.ident)
  PRiBAT_ASPC_cob$orig.ident <- as.factor(PRiBAT_ASPC_cob$orig.ident)
  saveRDS(PRiBAT_ASPC_cob, file = paste0(PRiBAT_ASPC_prefix,"Seurat.rds"))
  
  Idents(PRiBAT_ASPC_cob) <- PRiBAT_ASPC_cob$integrated_snn_res.0.1
  DimPlot(PRiBAT_ASPC_cob, reduction = "umap", label  = F, repel = T, cols = hughie_color, label.box = T, split.by = "orig.ident")
  Plot_Cell_compoistion(Seurat_Obj = PRiBAT_ASPC_cob, ColorUse1 = hughie_color, ColorUse2 = hughie_color, OutPrefix = PRiBAT_ASPC_prefix)
  TopMarkersGenes <- MarkersPlot(Seurat_Obj = PRiBAT_ASPC_cob, ColorUse = hughie_color, OutPrefix = PRiBAT_ASPC_prefix)
}

PRiBAT_ASPC_cob <- RenameIdents(PRiBAT_ASPC_cob, `0` = "aspc2", `1` = "aspc1", `2` = "aspc3")
PRiBAT_ASPC_cob$Myannotation_ASPC <- PRiBAT_ASPC_cob@active.ident
PRiBAT_ASPC_Ranks <- c("aspc1","aspc2","aspc3")
PRiBAT_ASPC_colors <- c("#0072B5FF","#E18727FF", "#6A3D9A")
names(PRiBAT_ASPC_colors) <- PRiBAT_ASPC_Ranks
Idents(PRiBAT_ASPC_cob) <- factor(Idents(PRiBAT_ASPC_cob), levels = PRiBAT_ASPC_Ranks)

PlotGenes <- c("Pdgfra","Dpp4","Pi16","Anxa3","Pdgfrb","Pparg","Cd36","Cidea","Ucp1")
VlnPlot(PRiBAT_ASPC_cob, features = PlotGenes, fill.by = "ident", cols = hughie_color, stack = T, flip=T, assay = "RNA")

(p1 <- VlnPlot(PRiBAT_ASPC_cob, features = PlotGenes, fill.by = "ident", cols = PRiBAT_ASPC_colors, stack = T, flip=T, assay = "RNA") + FontSize(x.text = 16, y.text = 16) + NoLegend())
(p2 <- VlnPlot(PRiBAT_ASPC_cob, features = PlotGenes, stack = T, flip=T, split.by = "orig.ident", cols = c("#729E4E","#6A3D9A"), assay = "RNA") + FontSize(x.text = 16, y.text = 16))
p1 + p2
ggsave(paste0(PRiBAT_ASPC_prefix, "Selected_markers_VlnPlot_", Sys.Date(), ".pdf"), width = 10, height = 8)

(p1 <- DimPlot(PRiBAT_ASPC_cob, reduction = "umap", group.by = c("ident"), label = F, repel = T, label.box =T, cols = PRiBAT_ASPC_colors, label.color = "white"))
(p2 <- DimPlot(PRiBAT_ASPC_cob, reduction = "umap", group.by = c("orig.ident"), label.box =T, cols = c("#729E4E","#6A3D9A")))
p1 + p2
ggsave(paste0(PRiBAT_ASPC_prefix, "UMAP-ByClusters_", Sys.Date(),".pdf"), width = 10, height = 5)

DimPlot(PRiBAT_ASPC_cob, reduction = "umap", label = T, split.by = "orig.ident", label.box =T, cols = PRiBAT_ASPC_colors, label.color = "white", repel = T,order = T)
ggsave(paste0(PRiBAT_ASPC_prefix,"UMAP-Byorig.ident_",Sys.Date(),".pdf"), width = 10, height = 5)

Plot_Cell_compoistion(Seurat_Obj = PRiBAT_ASPC_cob, OutPrefix = PRiBAT_ASPC_prefix, ColorUse1 = PRiBAT_ASPC_colors, ColorUse2 = c("#729E4E","#6A3D9A"))
TopMarkersGenes <- MarkersPlot(Seurat_Obj = PRiBAT_ASPC_cob, ColorUse = PRiBAT_ASPC_colors, OutPrefix = PRiBAT_ASPC_prefix,NMarkers = 20)
DEG_enrichment(Seruat_DEG_file = paste0(PRiBAT_ASPC_prefix, "Allmarkers_",Sys.Date(),".csv"), showCategoryNum = 20, filterLevel = 4, cellRanks = PRiBAT_ASPC_Ranks)

FeaturePlot(PRiBAT_ASPC_cob, features = c("Pdgfra","Dpp4","Pdgfrb","Cidea"), min.cutoff = "q10", order = T, pt.size=0.1, ncol = 2)
ggsave(paste0(PRiBAT_ASPC_prefix,"FeaturePlot_",Sys.Date(),".pdf"), width = 6, height = 6)
saveRDS(PRiBAT_ASPC_cob, file = paste0(PRiBAT_ASPC_prefix,"Seurat.rds"))

