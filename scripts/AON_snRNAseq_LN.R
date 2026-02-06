library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(scCustomize)
library(presto)
library(magrittr)
library(patchwork)
library(qs)
#library(Matrix)

AON.data <- Read10X_h5("../data/aon_10x/filtered_feature_bc_matrix.h5")

AON <- CreateSeuratObject(counts = AON.data, project = "AON3k", min.cells = 3, min.features = 200)
AON

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AON <- PercentageFeatureSet(AON, pattern = "^MT-", col.name = "percent.mt")

# Visualize QC metrics as a violin plot
VlnPlot(AON, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(AON, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AON, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

AON <- subset(AON, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

AON <- SCTransform(AON, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
AON <- RunPCA(AON, verbose = FALSE)
AON <- RunUMAP(AON, dims = 1:30, verbose = FALSE)

AON <- FindNeighbors(AON, dims = 1:30, verbose = FALSE)
AON <- FindClusters(AON, verbose = FALSE)
DimPlot(AON, label = FALSE)

VlnPlot(AON, features = c("Gad1", "Npy"),
        pt.size = 0.2, ncol = 2)

VlnPlot(AON, features = c( "Slc44a1", "Chat"),
        pt.size = 0.2, ncol = 2)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AON), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AON)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AON)

# Look at cluster IDs of the first 5 cells
head(Idents(AON), 5)

AON.markers <- FindAllMarkers(AON, only.pos = TRUE)
AON.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# saveRDS(AON, file = "/home/siddharth/AON_seq/AON_first_pass_v1.rds")

VlnPlot(AON, features = c("Slc17a7", "Camk2d"))
# you can plot raw counts as well
VlnPlot(AON, features = c("Slc17a7", "Camk2d"), slot = "counts", log = TRUE)

FeaturePlot(AON, features = c("Calb1", "Cck"), cols = c("lightgrey","purple"),pt.size = .5)

#Playing with scCustomize
candidate_genes <- read.csv("../data/aon_10x/Candidate_Genes.csv")

FeaturePlot_scCustom(seurat_object=AON,colors_use = viridis_dark_high, features=c("Ebf1","Ankk1"))

cluster_stats <- Cluster_Stats_All_Samples(seurat_object = AON)

cluster_stats

median_stats <- Median_Stats(seurat_object = AON)

all_markers_pct <- Add_Pct_Diff(marker_dataframe = AON.markers)

top_5 <- Extract_Top_Markers(marker_dataframe = AON.markers, num_genes = 5, rank_by = "avg_log2FC")


top_5_df <- Extract_Top_Markers(marker_dataframe = AON.markers, num_genes = 5, data_frame = TRUE,
                                rank_by = "avg_log2FC",make_unique=TRUE)

Iterate_FeaturePlot_scCustom(seurat_object = AON, features = candidate_genes$GeneName, single_pdf = T,
            file_path = "../output/",
            file_name = "Candidate_genes.pdf", features_per_page=4) & NoLegend()

Iterate_VlnPlot_scCustom(seurat_object = AON, features = candidate_genes$GeneName, single_pdf = T,
                             file_path = "../output/",
                             file_name = "Candidate_genes_medians.pdf",group.by="seurat_clusters", plot_median = TRUE) & NoLegend()


DotPlot_scCustom(seurat_object = AON, features = c("Slc17a6","Slc17a7" ,"Emx1", "Gad1","Calb1","Calb2","Cck","Vip","Npy","Pvalb", "Sst","Drd2", "Drd3",
                                                  "Gabbr1","Htr1d","Htr2a","Htr7" ),x_lab_rotate = TRUE,dot.scale= 10)

DotPlot_scCustom(seurat_object = AON, features = c("Slc17a6","Slc17a7", "Gad1","Calb1","Calb2","Cck","Vip","Npy","Pvalb", "Sst" ),x_lab_rotate = TRUE,dot.scale= 10)

ExcitatoryClusters <- subset(x = AON, subset = seurat_clusters == c("23", "22", "18","17","14","13"))
InhibitoryClusters <- subset(x = AON, subset = seurat_clusters == c("25", "21", "20","19","16","15","12","11","10","9","8","7","6","5","3","2","1","0"))


DotPlot_scCustom(seurat_object = InhibitoryClusters, features = c("Slc32a1","Gad1","Calb1","Calb2","Cck","Vip","Npy","Pvalb", "Sst","Camk2d","Dcx", "Drd2", "Drd3",
                                                   "Gabbr1","Htr1d","Htr2a","Htr7" ),x_lab_rotate = TRUE,dot.scale= 10)

DotPlot_scCustom(seurat_object = ExcitatoryClusters, features = c("Slc17a6","Slc17a7","Drd2", "Drd3","Gabbr1","Htr1d","Htr2a","Htr7","Fam81a","Ndst4","Nr4a2","Ociad2","Scg3","Vxn"),x_lab_rotate = TRUE,dot.scale= 10)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(AON, ident.1 = 2)
head(cluster2.markers, n = 5)


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(AON, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones

cluster0.markers <- FindMarkers(AON, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)




FeaturePlot(AON, features = c("Nmbr", "Ndst4"))
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

