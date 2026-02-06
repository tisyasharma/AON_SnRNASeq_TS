library(dplyr)
library(Seurat)
library(patchwork)
#library(Matrix)
AON.data <- Read10X_h5("../data/aon_10x/filtered_feature_bc_matrix.h5")

AON <- CreateSeuratObject(counts = AON.data, project = "AON3k", min.cells = 3, min.features = 200)
AON

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AON[["percent.mt"]] <- PercentageFeatureSet(AON, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(AON, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(AON, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AON, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

AON <- subset(AON, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

AON <- NormalizeData(AON, normalization.method = "LogNormalize", scale.factor = 10000)

AON<- NormalizeData(AON)
AON <- FindVariableFeatures(AON, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AON), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AON)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(AON)
AON <- ScaleData(AON, features = all.genes)


AON <- RunPCA(AON, features = VariableFeatures(object = AON))

# Examine and visualize PCA results a few different ways
print(AON[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(AON, dims = 1:2, reduction = "pca")

DimPlot(AON, reduction = "pca") + NoLegend()

DimHeatmap(AON, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(AON, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(AON)

AON <- FindNeighbors(AON, dims = 1:10)
AON <- FindClusters(AON, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(AON), 5)

AON <- RunUMAP(AON, dims = 1:10)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(AON, reduction = "umap")

saveRDS(AON, file = "../data/aon_10x/AON_first_pass_v1.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(AON, ident.1 = 2)
head(cluster2.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(AON, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
AON.markers <- FindAllMarkers(AON, only.pos = TRUE)
AON.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(AON, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


VlnPlot(AON, features = c("Slc17a7", "Camk2d"))
# you can plot raw counts as well
VlnPlot(AON, features = c("Slc17a7", "Camk2d"), slot = "counts", log = TRUE)

FeaturePlot(AON, features = c("Slc17a7", "Slc32a1"), cols = c("lightgrey","forestgreen"),pt.size = .5)

FeaturePlot(AON, features = c("Slc17a7", "Gria2"))
FeaturePlot(AON, features = c("Slc17a7", "Grin2b"))
FeaturePlot(AON, features = c("Rbfox3", "Nts","Uts2b","Nucb2"))
FeaturePlot(AON, features = c("Rbfox3", "Calca","Adcyap1","Calcb"))
FeaturePlot(AON, features = c("Rbfox3", "Avp","Pdyn","Apln"))
FeaturePlot(AON, features = c("Rbfox3", "Pomc","Npvf","Grp"))
FeaturePlot(AON, features = c("Rbfox3", "Gcg","Pnoc","Cart"))
FeaturePlot(AON, features = c("Rbfox3", "Ucn","Crh","Penk"))
FeaturePlot(AON, features = c("Rbfox3", "Tac1","Nucb2","Penk"))        #neuropeptides in PVT expressed sourced from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873951/


FeaturePlot(AON, features = c("Slc17a7", "Slc32a1"), cols = c("lightgrey","blue"),pt.size = .2)
FeaturePlot(AON, features = c("Camk2d", "Gad2"), cols = c("lightgrey","blue"),pt.size = .2)

FeaturePlot(AON, features = c("Reln", "Slc17a7"), cols = c("lightgrey","blue"),pt.size = .2)

FeaturePlot(AON, features = c("Syt1", "Slc17a7"), cols = c("lightgrey","blue"),pt.size = .2)

  FeaturePlot(AON,features = c("Sst", "Pthlh", "Scg3", "Nxph3" ))
FeaturePlot(AON,features = c("Sst", "Pvalb", "Calr", "Calb1" ))




## TO DO

IEG.genes = c("Fos","Npas4","Arc","Jun","Egr1")

percent.IEG = colSums(as.array(expm1(AON@data[IEG.genes, ])))/colSums(as.array(expm1(AON@data)))
AON_v1 = AddMetaData(AON, percent.IEG, "percent.IEG")

VlnPlot(OB, c("percent.RBC", "percent.apoptotic","percent.IEG"))




n.cells <- 1000 # random sampling from 1000 cells
AON_1 <- subset(pbmc, cells.use = sample(pbmc@cell.names, size = n.cells))
length(AON_1@cell.names)
OB_RS_1 = FindClusters(OB_RS_1, genes.use = OB_RS_1@var.genes, dims.use = 1:40, 
                       k.param = 30, save.SNN = TRUE,
                       algorithm = 3, resolution = 1, force.recalc = T)
OB_RS_1 = RunTSNE(OB_RS_1, reduction.use = "pca", dims.use = 1:40, max_iter = 10000, perplexity = 20, theta = 0.25, verbose = TRUE)

TSNEPlot(object = OB_RS_1, do.label = T)