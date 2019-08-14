### TITLE : Integrated analysis with Seurat
### AUTHOR : Javier Perales-Patón - javier.perales@bioquant.uni-heidelberg.de
### DESCRIPTION : Integrated analysis with Seurat to for two scRNA-seq samples from two different conditions


suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(GSEABase))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(clustree))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(optparse))
# Get some functions for Seurat analysis
source("./src/seurat_fx.R")

### 2 Get input parameters #####
#--- Input variables
option_list = list(
  make_option(c("--CASE_INPUT"), action="store", default="./data/sc/disease/filtered_feature_bc_matrix/", type='character',
              help="cellranger count folder that contains the output matrices"),
  make_option(c("--CASE_SNAME"), action="store", default="disease", type='character',
              help="Sample name"),
  make_option(c("--CONTROL_INPUT"), action="store", default="./data/sc/control/filtered_feature_bc_matrix/", type='character',
              help="cellranger count folder that contains the output matrices"),
    make_option(c("--CONTROL_SNAME"), action="store", default="control", type='character',
              help="Sample name"),
  make_option(c("--OUTDIR"), action="store", default="./results/Integrated_disease_control/", type='character',
              help="Output directory"),
  make_option(c("--NPC"), action="store", default=50, type='numeric',
              help="number of principal components to calculate."),
  make_option(c("--NPC_ANCHOR"), action="store", default=20, type='numeric',
              help="number of principal components to consider for the anchoring."),
  make_option(c("--NPC_CLUSTERING"), action="store", default=20, type='numeric',
              help="Number of principal components to consider for the cell clustering."),
  make_option(c("--RES"), action="store", default=0.3, type='numeric',
              help="Resolution for cell clustering.")
)

# Parse the parameters
opt = parse_args(OptionParser(option_list=option_list))

# Cat the input parameters
cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}


### 3 XXXXXXX ####
# Read samples and estimate perc mitocondrial genes
CASE <- getSeuratObject(path = CASE_INPUT, project_name = CASE_SNAME, mt.pattern = "^MT-", min.cells = 3, min.features = 200)
CASE$stim <- CASE_SNAME
CASE <- subset(CASE, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 80)
CASE <- NormalizeData(CASE, verbose = FALSE)
CASE <- FindVariableFeatures(CASE, selection.method = "vst", nfeatures = 2000)
  
CTRL <- getSeuratObject(path = CONTROL_INPUT, project_name = CONTROL_SNAME, mt.pattern = "^MT-", min.cells = 3, min.features = 200)
CTRL$stim <- CONTROL_SNAME
CTRL <- subset(CTRL, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 80)
CTRL <- NormalizeData(CTRL, verbose = FALSE)
CTRL <- FindVariableFeatures(CTRL, selection.method = "vst", nfeatures = 2000)

### 4 Perform integration #####
anchors <- FindIntegrationAnchors(object.list = list(CTRL, CASE), dims = 1:NPC_ANCHOR)
S <- IntegrateData(anchorset = anchors, dims = 1:20)

# Get some RAM space: we don't need this anymore
rm(CASE, CTRL)

# Create folder to store data
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

### 5 Perform an integrated analysis ####
DefaultAssay(S) <- "integrated"

# Run the standard workflow for visualization and clustering
S <- ScaleData(S, verbose = FALSE)
S <- RunPCA(S, npcs = NPC, verbose = FALSE)

png(paste0(OUTDIR,"/elbow_pca.png"),width = 600*3,height = 800*3, res=280)
ElbowPlot(S,ndims = NPC) + geom_vline(xintercept = NPC_CLUSTERING, col="red")
dev.off()

# t-SNE and Clustering
S <- FindNeighbors(S, reduction = "pca", dims = 1:NPC_CLUSTERING)
S <- FindClusters(S, resolution = seq(from=0.1, to=1.5, by=0.1))

# clustree(S, prefix = "RNA_snn_res.")
saveClusterTree(S, OUTDIR, prefix="integrated_snn_res.")

S <- FindClusters(S, resolution = RES)

# Run UMAP
S <- RunUMAP(S, reduction = "pca", dims = 1:NPC_CLUSTERING)

png(paste0(OUTDIR,"/umap.png"),width = 800*3, height = 700*3, res=280)
DimPlot(S, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 14)
dev.off()

# Visualization
p1 <- DimPlot(S, reduction = "umap", label = TRUE, label.size = 14)
p2 <- DimPlot(S, reduction = "umap", group.by = "stim")

png(paste0(OUTDIR,"/umap_integrated.png"),width = 1600*3, height = 700*3, res=280)
plot_grid(p1, p2)
dev.off()

### 6 Switch to RNA ####
DefaultAssay(S) <- "RNA"
# cl1.markers <- FindConservedMarkers(S, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
# head(cl1.markers)

MARKERS.OUTDIR <- paste0(OUTDIR,"/markers")
if(!dir.exists(MARKERS.OUTDIR)) dir.create(MARKERS.OUTDIR);

markers <- setNames(vector("list",length=length(levels(S))), levels(S))
for(idx in names(markers)) {
  markers.idx <- FindConservedMarkers(S,ident.1 = idx,ident.2 = setdiff(levels(S), idx),grouping.var = "stim",only.pos=T)
  cols_names <- colnames(markers.idx)
  
  # Add two extra cols
  markers.idx$cluster <- idx
  markers.idx$gene <- rownames(markers.idx)
  
  write.table(markers.idx[,c("cluster", "gene", cols_names)],
              file = paste0(MARKERS.OUTDIR,"/cluster",idx,".tsv"),
              sep="\t",col.names = TRUE, row.names = FALSE, quote=FALSE
  )
}


# S <- NormalizeData(S, verbose = FALSE) # I am not sure if I have to normalize again
S <- Seurat_scaledata(S)

### 7 Check markers
png(paste0(OUTDIR,"/heatmap_chen.png"),width = 1500*3,height = 1000*3, res=280)
DoHeatmap2(SeuratObject = S, GSC = getGmt("./REFERENCES/Chen_etal2019/Chen_JASN.gmt"),
           assay = "RNA", res=NULL, show_hr = FALSE)
dev.off()

png(paste0(OUTDIR,"/heatmap_Lake.png"),width = 1500*3,height = 1000*3, res=280)
DoHeatmap2(SeuratObject = S, GSC = getGmt("./REFERENCES/Lake_etal2019/markers.gmt"),
           assay = "RNA", res=NULL)
dev.off()

### 8 Diff density of populations ######
cl.cnt<-tapply(S$seurat_clusters,S@meta.data$stim,table)
# sapply(names(cl.cnt), function(idx) sum(S$stim==idx))
cl.prop <- cl.cnt
for(idx in names(cl.prop)) {
  cl.prop[[idx]] <- (cl.prop[[idx]]/sum(S$stim==idx)*100)
}

# DATA <- data.frame(matrix(0, ncol=3, nrow=length(levels(S))*2, dimnames=list(1:(length(levels(S))*2),c("Cluster","Sample","Percentage"))))
# DATA$Cluster <- rep(levels(S),2) 
# DATA$Sample <- unlist(sapply(names(cl.prop), function(z) rep(z, length=length(levels(S))), simplify = FALSE))
# stopifnot(all(DATA$Cluster== unlist(sapply(names(cl.prop), function(z) names(cl.prop[[z]]), simplify = FALSE))))
# DATA$Percentage <- unlist(cl.prop)

DATA <- data.frame(cluster=levels(S),
                   S1=as.numeric(cl.prop[[1]]),
                   S2=as.numeric(cl.prop[[2]]))
# colnames(DATA)[-1] <- names(cl.prop)
# DATA <- rev(DATA)

# Source: https://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library(grid)
g.mid<-ggplot(DATA,aes(x=1,y=cluster))+geom_text(aes(label=cluster), size=22)+
  geom_segment(aes(x=0.94,xend=0.96,yend=cluster))+
  geom_segment(aes(x=1.04,xend=1.065,yend=cluster))+
  ggtitle("Cluster")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(plot.title = element_text(size=28),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA, size=28),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))

g1 <- ggplot(data = DATA, aes(x = cluster, y = S1, fill=cluster)) +
  geom_bar(stat = "identity") + ggtitle(names(cl.cnt)[1]) +
  theme(plot.title = element_text(size=28),
        legend.position = "none",
        axis.text.x = element_text(size=28),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

g2 <- ggplot(data = DATA, aes(x = cluster, y = S2, fill=cluster)) +xlab(NULL)+
  geom_bar(stat = "identity") + ggtitle(names(cl.prop)[2]) +
  theme(plot.title = element_text(size=28),
        legend.position = "none",
        axis.text.x = element_text(size=28),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  coord_flip()

library(gridExtra)
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

png(paste0(OUTDIR,"/cluster_density.png"),width = 1000*3, height = 1000*3, res=280)
grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9),
             bottom = textGrob("Percentage of cells within-sample (%)",
                               gp=gpar(fontsize=32,font=2)))
dev.off()

### 9 Cell phase ####
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

S <- CellCycleScoring(S, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

png(paste0(OUTDIR,"/umap_cellcycle.png"),width = 800*3, height = 700*3, res=280)
DimPlot(S, reduction = "umap", group.by = "Phase")
dev.off()

### 10 Mt ###
# Add mitocondrial genes
S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")
d1 <- FeaturePlot(S, reduction = "umap", features = "percent.mt")
d2 <- FeaturePlot(S, reduction = "umap", features = "nCount_RNA")

png(paste0(OUTDIR,"/umap_diagnostic.png"),width = 1600*3, height = 700*3, res=280)
plot_grid(d1, d2)
dev.off()

### 11 Add disease+ : the marker of cell sorting ####
png(paste0(OUTDIR,"/umap_cd24marker.png"),width = 1600*3, height = 700*3, res=280)
FeaturePlot(S, reduction = "umap", features = "disease",split.by = "stim", label = TRUE, label.size = 14)
dev.off()

# ### 14 Hypothesis-driven tests ####
# png(paste0(OUTDIR,"/markers_immunecell.png"),width = 800*3,height = 700*3,res = 280)
# FeaturePlot(S, features = "Ptprc", label = TRUE, label.size = 14)
# # + ggtitle("immune cell marker: Ptprc")
# dev.off()
