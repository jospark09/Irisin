# 1. Setting ----
# work on file name with complementary base and reverse its orders

library(Matrix)
library(stringr)
library(Seurat)
library(purrr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tribe)
library(disk.frame)
library(patchwork)
library(stringi)
library(clusterProfiler)
library(xlsx)
library(BiocManager)
library(remotes)
library(devtools)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)
library(gplots)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)
theme_set(theme_cowplot())


# Color setting
my_color <- c('yellow', 'red')

# dotplot scale by size
dot <- function(obj, gene, title){
  DotPlot(obj, features = gene, dot.scale = 10, cols = my_color, scale.by = "size") +
    RotatedAxis() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5))  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}


# dotplot scale by radius
dot01 <- function(obj, gene, title){
  DotPlot(obj, features = gene, dot.scale = 10, cols = my_color, scale.by = "radius") +
    RotatedAxis() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5))  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}

barcodes.list <- list("CTATTAAG", 
              "AAGGCTAT",
              "GAGCCTTA",
              "TTATGCGA",
              "GGAGGTAA",
              "CATAACTG",
              "AGTAAAGG",
              "TCCGTCTC",
              "AGCTTTCT",
              "AAGAGCGT",
              "AGAATGCG")


# reverse the barcodes
barcodes_rv <- lapply(barcodes.list, stri_reverse)
unlist(barcodes.list)
unlist(barcodes_rv)
  
# converting function
convert <- function(x){
  chartr("ATGC","TACG",x)
}


# switch each letter to the mapped complementary base
converted_barcodes <- lapply(barcodes_rv, convert)
unlist(converted_barcodes)


# 2. Read-in files ----
library(Matrix)
read.in <- function(x){read.table(x, as.is = T, header = F)[, 1]}

# first batch
h10_01_a <- readMM("a_h10_01_Eunhee_Drug_Treat-CTTAATAG/Eunhee_Drug_Treat-CTTAATAG.mtx")
h10_01_a.cn <- read.in("a_h10_01_Eunhee_Drug_Treat-CTTAATAG/Eunhee_Drug_Treat-CTTAATAG.mtx.colnames")
h10_01_a.rn <- read.in("a_h10_01_Eunhee_Drug_Treat-CTTAATAG/Eunhee_Drug_Treat-CTTAATAG.mtx.rownames")

h10_03_a <- readMM("a_h10_03_Eunhee_Drug_Treat-ATAGCCTT/Eunhee_Drug_Treat-ATAGCCTT.mtx")
h10_03_a.cn <- read.in("a_h10_03_Eunhee_Drug_Treat-ATAGCCTT/Eunhee_Drug_Treat-ATAGCCTT.mtx.colnames")
h10_03_a.rn <- read.in("a_h10_03_Eunhee_Drug_Treat-ATAGCCTT/Eunhee_Drug_Treat-ATAGCCTT.mtx.rownames")

a5_01_a <- readMM("a_a5_01_Eunhee_Drug_Treat-TAAGGCTC/Eunhee_Drug_Treat-TAAGGCTC.mtx")
a5_01_a.cn <- read.in("a_a5_01_Eunhee_Drug_Treat-TAAGGCTC/Eunhee_Drug_Treat-TAAGGCTC.mtx.colnames")
a5_01_a.rn <- read.in("a_a5_01_Eunhee_Drug_Treat-TAAGGCTC/Eunhee_Drug_Treat-TAAGGCTC.mtx.rownames")

a5_02_a <- readMM("a_a5_02_Eunhee_Drug_Treat-TCGCATAA/Eunhee_Drug_Treat-TCGCATAA.mtx")
a5_02_a.cn <- read.in("a_a5_02_Eunhee_Drug_Treat-TCGCATAA/Eunhee_Drug_Treat-TCGCATAA.mtx.colnames")
a5_02_a.rn <- read.in("a_a5_02_Eunhee_Drug_Treat-TCGCATAA/Eunhee_Drug_Treat-TCGCATAA.mtx.rownames")

a5_03_a <- readMM("a_a5_03_Eunhee_Drug_Treat-TTACCTCC/Eunhee_Drug_Treat-TTACCTCC.mtx")
a5_03_a.cn <- read.in("a_a5_03_Eunhee_Drug_Treat-TTACCTCC/Eunhee_Drug_Treat-TTACCTCC.mtx.colnames")
a5_03_a.rn <- read.in("a_a5_03_Eunhee_Drug_Treat-TTACCTCC/Eunhee_Drug_Treat-TTACCTCC.mtx.rownames")


# second batch
h10_01_b <- readMM("b_h10_01_Eunhee_Drug_Treat-AGAAAGCT/Eunhee_Drug_Treat-AGAAAGCT.mtx")
h10_01_b.cn <- read.in("b_h10_01_Eunhee_Drug_Treat-AGAAAGCT/Eunhee_Drug_Treat-AGAAAGCT.mtx.colnames")
h10_01_b.rn <- read.in("b_h10_01_Eunhee_Drug_Treat-AGAAAGCT/Eunhee_Drug_Treat-AGAAAGCT.mtx.rownames")

h10_02_b <- readMM("b_h10_02_Eunhee_Drug_Treat-ACGCTCTT/Eunhee_Drug_Treat-ACGCTCTT.mtx")
h10_02_b.cn <- read.in("b_h10_02_Eunhee_Drug_Treat-ACGCTCTT/Eunhee_Drug_Treat-ACGCTCTT.mtx.colnames")
h10_02_b.rn <- read.in("b_h10_02_Eunhee_Drug_Treat-ACGCTCTT/Eunhee_Drug_Treat-ACGCTCTT.mtx.rownames")

h10_03_b <- readMM("b_h10_03_Eunhee_Drug_Treat-CGCATTCT/Eunhee_Drug_Treat-CGCATTCT.mtx")
h10_03_b.cn <- read.in("b_h10_03_Eunhee_Drug_Treat-CGCATTCT/Eunhee_Drug_Treat-CGCATTCT.mtx.colnames")
h10_03_b.rn <- read.in("b_h10_03_Eunhee_Drug_Treat-CGCATTCT/Eunhee_Drug_Treat-CGCATTCT.mtx.rownames")

a5_01_b <- readMM("b_a5_01_Eunhee_Drug_Treat-CAGTTATG/Eunhee_Drug_Treat-CAGTTATG.mtx")
a5_01_b.cn <- read.in("b_a5_01_Eunhee_Drug_Treat-CAGTTATG/Eunhee_Drug_Treat-CAGTTATG.mtx.colnames")
a5_01_b.rn <- read.in("b_a5_01_Eunhee_Drug_Treat-CAGTTATG/Eunhee_Drug_Treat-CAGTTATG.mtx.rownames")

a5_02_b <- readMM("b_a5_02_Eunhee_Drug_Treat-CCTTTACT/Eunhee_Drug_Treat-CCTTTACT.mtx")
a5_02_b.cn <- read.in("b_a5_02_Eunhee_Drug_Treat-CCTTTACT/Eunhee_Drug_Treat-CCTTTACT.mtx.colnames")
a5_02_b.rn <- read.in("b_a5_02_Eunhee_Drug_Treat-CCTTTACT/Eunhee_Drug_Treat-CCTTTACT.mtx.rownames")

a5_03_b <- readMM("b_a5_03_Eunhee_Drug_Treat-GAGACGGA/Eunhee_Drug_Treat-GAGACGGA.mtx")
a5_03_b.cn <- read.in("b_a5_03_Eunhee_Drug_Treat-GAGACGGA/Eunhee_Drug_Treat-GAGACGGA.mtx.colnames")
a5_03_b.rn <- read.in("b_a5_03_Eunhee_Drug_Treat-GAGACGGA/Eunhee_Drug_Treat-GAGACGGA.mtx.rownames")


# make a list of datasets
data.list <- list(h10_01_a,
                  h10_03_a,
                  h10_01_b,
                  h10_02_b,
                  h10_03_b,
                  a5_01_a,
                  a5_02_a,
                  a5_03_a,
                  a5_01_b,
                  a5_02_b,
                  a5_03_b)

names(data.list) <- c("h10_01_a",
                      "h10_03_a",
                      "h10_01_b",
                      "h10_02_b",
                      "h10_03_b",
                      "a5_01_a",
                      "a5_02_a",
                      "a5_03_a",
                      "a5_01_b",
                      "a5_02_b",
                      "a5_03_b")

unlist(lapply(data.list, nrow))
unlist(lapply(data.list, ncol))


# Gene name annotation ----
# before doing BioMart, check all rownames to see they are all identical

# check if all rn vectors are identical
rn.list <- list(h10_01_a.rn,
                h10_03_a.rn,
                h10_01_b.rn,
                h10_02_b.rn,
                h10_03_b.rn,
                a5_01_a.rn,
                a5_02_a.rn,
                a5_03_a.rn,
                a5_01_b.rn,
                a5_02_b.rn,
                a5_03_b.rn)


for (i in 1:(length(rn.list)-1)){
  print(identical(rn.list[[i]], rn.list[[i+1]]))
}
# All TRUE

# getBM
# Optimize memory space before processing and split the object by sample for iterative process
options(future.globals.maxSize = 4000 * 1024^2)
library(biomaRt)
# Retrieve BioMart
ensembl <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "www.ensembl.org")

# hgnc symbols matching with rownames
bm <- biomaRt::getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = h10_01_a.rn,
                     mart= ensembl)

# BioMart cleaning
# 01. empty symbol
ids <- which(nchar(bm$hgnc_symbol)==0)
length(ids)     # 18059
bm <- bm[-ids,]

# 02. duplicated
ids.dup <- which(duplicated(bm$hgnc_symbol))
length(ids.dup)     # 17
bm <- bm[-ids.dup,]

# 03. final indices of ensembl ids (good.ids will be used for all rownames)
good.ids <- match(bm$ensembl_gene_id, h10_01_a.rn)

# check
identical(bm$ensembl_gene_id, h10_01_a.rn[good.ids])


# _ Assembly ----
## 01. Select only rows that are good.ids in each dataset
## 02. Assemble rownames: rownames are hgnc column of the cleaned bm
## 03. Assemble colnames


# Assembly function
cleaning <- function(x){
  x <- x[good.ids, ]
}


good.data.list <- list()
good.data.list <- lapply(data.list, cleaning)  # good.data.list is the list of data with clean row names (no NAs, no duplicates in hgnc symbol)

# rownames
for(i in 1:length(good.data.list)){
  rownames(good.data.list[[i]]) <- bm$hgnc_symbol
}

# colnames
cn.list <- list(h10_01_a.cn,
                h10_03_a.cn,
                h10_01_b.cn,
                h10_02_b.cn,
                h10_03_b.cn,
                a5_01_a.cn,
                a5_02_a.cn,
                a5_03_a.cn,
                a5_01_b.cn,
                a5_02_b.cn,
                a5_03_b.cn)


for(i in 1:length(good.data.list)){
  colnames(good.data.list[[i]]) <- cn.list[[i]]
}


# Check: no change in column numbers after cleaning. cleaning was only for rows.
identical(unlist(lapply(good.data.list, ncol)), unlist(lapply(data.list, ncol)))

# 3. Seurat object ----
seurat.obj <- function(x, n){CreateSeuratObject(counts = x, 
                                             project = "Irisin", 
                                             min.cells = 3, min.features = n)}






# make Seurat objects
library(Seurat)
seurat.list <- list()
seurat.list <- lapply(good.data.list, seurat.obj, n=200)    # selected for dowbstream analysis

# cell counts 
unlist(lapply(seurat.list, ncol))

# working with min.features = 200
# add metadata
for(i in 1:length(seurat.list)){
  seurat.list[[i]]$sample <- names(seurat.list[i])
  seurat.list[[i]]$tech <- "inDrop"
  seurat.list[[i]]$stim <- "AD"
  seurat.list[[i]]$batch <- "A"
  seurat.list[[i]]$tx <- "Irisin + Anti.Integrin"
  seurat.list[[i]]$line <- "H10"
}


# Edit batch
for(i in c(3:5, 9:11)){
  seurat.list[[i]]$batch <- "B"
}

# Edit tx (_01: None, _02: Tx1, 03_IGG)
for(i in c(1,3, 6, 9)){
  seurat.list[[i]]$tx <- "None"
}

for(i in c(4,7, 10)){
  seurat.list[[i]]$tx <- "Irisin"
}


# Edit line
for(i in c(6:11)){
  seurat.list[[i]]$line <- "A5"
}



for(i in 1:11){
  print(table(seurat.list[[i]]$line))
  print(table(seurat.list[[i]]$batch))
  print(table(seurat.list[[i]]$tx))
}


# Add meta mito content
for(i in 1:length(seurat.list)){
  seurat.list[[i]]$percent.mt <- PercentageFeatureSet(seurat.list[[i]], pattern = "^MT-")
}

for(i in 1:length(seurat.list)){
  seurat.list[[i]] <- subset(x = seurat.list[[i]],
                             subset= (nFeature_RNA >= 200) &
                               (percent.mt <= 20))
}

# check all metadata items
names(seurat.list[[1]]@meta.data)


# VlnPlot Function setting
violin <- function(seurat){
  VlnPlot(seurat,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3, 
          pt.size = 0)
}


pdf(file = "figures/violin.pdf")
  violin(seurat.list[[1]])
  violin(seurat.list[[2]])
  violin(seurat.list[[3]])
  violin(seurat.list[[4]])
  violin(seurat.list[[5]])
  violin(seurat.list[[6]])
  violin(seurat.list[[7]])
  violin(seurat.list[[8]])
  violin(seurat.list[[9]])
  violin(seurat.list[[10]])
  violin(seurat.list[[11]])
dev.off()





# # _ QC visualization ----
# # Cell Counts by Sample
# library(ggplot2)
# library(tibble)
# metadata %>% 
#   ggplot(aes(x=sample, fill=sample)) + 
#   geom_bar() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells")
# 
# 
# # Number of UMIs/transcripts per cell
# metadata %>% 
#   ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   scale_x_log10() + 
#   theme_classic() +
#   ylab("Cell density") +
#   geom_vline(xintercept = 300)
# 
# 
# # Detected genes per cell: distribution
# metadata %>% 
#   ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
#   geom_density(alpha = 0.2) + 
#   theme_classic() +
#   scale_x_log10() + 
#   geom_vline(xintercept = c(150, 200, 300)) +
#   ggtitle("Unique Genes per Cell") +
#   theme(plot.title = element_text(hjust=0.5, face="bold"))
# 
# 
# # Detected genes per cell: distribution via boxplot
# metadata %>% 
#   ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
#   geom_boxplot() + 
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("nFeature_RNAs per Cell")
# 
# # Number of genes per UMI
# # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
# metadata %>% 
#   ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
#   geom_point() + 
#   scale_colour_gradient(low = "gray90", high = "black") +
#   stat_smooth(method=lm) +
#   scale_x_log10() + 
#   scale_y_log10() + 
#   theme_classic() +
#   geom_vline(xintercept = 300) +
#   geom_hline(yintercept = 200) +
#   facet_wrap(~sample)
# 
# 
# # Complexity of the gene expression: Detected genes per UMI
# metadata %>%
#   ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
#   geom_density(alpha = 0.2) +
#   theme_classic() +
#   geom_vline(xintercept = 0.8)

# Integration ----

for (i in 1:11) {
  seurat.list[[i]] <- SCTransform(seurat.list[[i]], vars.to.regress = c("percent.mt"))
}


features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)

seurat.list <- PrepSCTIntegration(object.list = seurat.list, 
                                 anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seurat.list,  normalization.method = "SCT", anchor.features = features)

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors,  normalization.method = "SCT")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering

combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

library(harmony)
combined.harmony <- combined %>% RunHarmony("sample", assay.use = "SCT", plot_convergence = T)

tri.harmony_embeddings <- Embeddings(combined.harmony, 'harmony')
tri.harmony_embeddings[1:5, 1:5]

DimPlot(object = combined.harmony, reduction = "harmony", pt.size = 2, split.by = "sample") + NoLegend()

DefaultAssay(combined.harmony) <- "SCT"
combined.harmony <- RunUMAP(combined.harmony, reduction = "harmony", dims = 1:30)
combined.harmony <- FindNeighbors(combined.harmony, reduction = "harmony", dims = 1:30)
combined.harmony <- FindClusters(combined.harmony, resolution = c(0.4, 0.8, 1.2, 1.6))

names(combined.harmony@meta.data)


Tri.res.list <- list(names(combined.harmony@meta.data[13]),
                     names(combined.harmony@meta.data[14]),
                     names(combined.harmony@meta.data[15]),
                     names(combined.harmony@meta.data[16]))

Tri.res.list

library(stringr)
names(Tri.res.list) <- str_sub(names(combined.harmony@meta.data[13:16]), start = 1L, end = -1L)
View(Tri.res.list)

# :: View clusters
# make sure the default assay is integrated
# DefaultAssay(seurat.integrated.combined) <- "SCT"
# DefaultAssay(seurat.integrated.combined) <- "integrated"

combined.harmony$SCT_snn_res.0.4 <- factor(combined.harmony$SCT_snn_res.0.4,
                                               levels = 0:(length(unique(combined.harmony$SCT_snn_res.0.4))-1))

combined.harmony$SCT_snn_res.0.8 <- factor(combined.harmony$SCT_snn_res.0.8,
                                               levels = 0:(length(unique(combined.harmony$SCT_snn_res.0.8))-1))

combined.harmony$SCT_snn_res.1.2 <- factor(combined.harmony$SCT_snn_res.1.2,
                                               levels = 0:(length(unique(combined.harmony$SCT_snn_res.1.2))-1))

combined.harmony$SCT_snn_res.1.6 <- factor(combined.harmony$SCT_snn_res.1.6,
                                               levels = 0:(length(unique(combined.harmony$SCT_snn_res.1.6))-1))

Idents(combined.harmony) <- "SCT_snn_res.0.4"

# combined.harmony$stim <- factor(combined.harmony$stim, levels = c("CTRL",  "AD"))
table(combined.harmony@meta.data[["sample"]])
combined.harmony$sample <- factor(combined.harmony$sample, 
                                   levels = c("a5_01_a", "a5_02_a", "a5_03_a",
                                              "a5_01_b", "a5_02_b", "a5_03_b",
                                              "h10_01_a", "h10_03_a",
                                              "h10_01_b", "h10_02_b", "h10_03_b"
                                              ))

DimPlot(combined.harmony,  label = TRUE, pt.size = 0.8, label.size = 7) + NoLegend()
DimPlot(combined.harmony, label = TRUE, pt.size = 0.8, label.size = 7, split.by = "sample") + NoLegend()
DimPlot(combined.harmony,  label = TRUE, pt.size = 0.8, label.size = 7, split.by = "each") + NoLegend()

table(combined.harmony@active.ident, combined.harmony@meta.data$tri)
saveRDS(combined.harmony, file = "objects/combined.harmony.mt20.RNA400.harmony.rds")
combined.harmony <- readRDS(file = "objects/combined.harmony.mt20.RNA400.harmony.rds")


DefaultAssay(combined.harmony) <- "SCT"
Trimarkers <- c( "GFAP", "VIM", "SOX2",  "AQP4", "ID4", "AGT",  "WSB1", "KCNQ1OT1", "KCNJ6", "SYT1",  "TOP2A", "MKI67")
dot01(combined.harmony, gene = Trimarkers, "Percent expression      Average gene expression")
dot(combined.harmony , gene = Trimarkers, "DotPlot (2.0)_by size")


DefaultAssay(combined.harmony) <- "RNA"
combined.harmony.Allmarkers <- FindAllMarkers(combined.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.harmony.Allmarkers.top5 <- combined.harmony.Allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

View(combined.harmony.Allmarkers.top5)
View(combined.harmony.Allmarkers)

library(xlsx)
write.xlsx(combined.harmony.Allmarkers, 
           file = "results/Find.all.markers.11.11.22.xlsx")

table(combined.harmony$stim)
res04.t <- t(as.data.frame(table(combined.harmony$integrated_snn_res.0.5)))
View(res04.t)
FeaturePlot(combined.harmony, pt.size = 2.0, features = c("CSF1R"))

table(combined.harmony@active.ident, combined.harmony@meta.data$tx)
composition <- data.frame(table(combined.harmony@active.ident, combined.harmony@meta.data$tx))
write.xlsx(composition, 
           file = "results/composition.10.14.22.xlsx")

res05.t <- t(as.data.frame(table(combined.harmony@active.ident, combined.harmony@meta.data$stim)))
View(res05.t)
table(combined.harmony@active.ident, combined.harmony@meta.data$each)

table(combined.harmony@meta.data[["sample"]])
Idents(combined.harmony) <- "sample"
Idents(seurat_integrated.3D.ast.a5.composition) <- "SCT_snn_res.0.4"
seurat_integrated.3D.ast.a5.composition <- subset(combined.harmony, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b", "a5_03_a", "a5_03_b"))

table(seurat_integrated.3D.ast.a5.composition@active.ident, seurat_integrated.3D.ast.a5.composition@meta.data$tx)
composition <- table(seurat_integrated.3D.ast.a5.composition@active.ident, seurat_integrated.3D.ast.a5.composition@meta.data$tx)
write.xlsx(composition, 
           file = "results/composition.a5.astrocyte.11.11.22.xlsx")



# 2022.11.11 3D name & DE ----
seurat_integrated.selected.01.combined.name <- RenameIdents(object = combined.harmony, 
                                                         "0" = "Astrocyte",
                                                         "1" = "Astrocyte",
                                                         "2" =  "RPS",
                                                         "3" =  "Neuron",
                                                         "4" = "Neuron",
                                                         "5" =  "Astrocyte",
                                                         "6" = "Proliferation",
                                                         "7" = "Astrocyte")


seurat_integrated.selected.01.combined.name.02 <- RenameIdents(object = combined.harmony, 
                                                            "0" = "Astrocyte1",
                                                            "1" = "Astrocyte2",
                                                            "2" =  "RPS",
                                                            "3" =  "Neuron1",
                                                            "4" = "Neuron2",
                                                            "5" =  "Astrocyte3",
                                                            "6" = "Proliferation",
                                                            "7" = "Astrocyte4")
seurat_integrated.selected.01.combined.name.02.subset <- subset(seurat_integrated.selected.01.combined.name.02, idents = c("Neuron1", "Neuron2", "Astrocyte1", "Astrocyte2", "Astrocyte3", "Astrocyte4"))

DimPlot(seurat_integrated.selected.01.combined.name.02.subset, label = TRUE, pt.size = 3, label.size = 7) + NoLegend()

DefaultAssay(seurat_integrated.selected.01.combined.name.02.subset) <- "RNA"
seurat_integrated.selected.01.combined.name.02.subset.Allmarkers <- FindAllMarkers(seurat_integrated.selected.01.combined.name.02.subset, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
seurat_integrated.selected.01.combined.name.02.subset.Allmarkers.top5 <- seurat_integrated.selected.01.combined.name.02.subset.Allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 40, order_by = avg_log2FC)


# $$ Gene.name.label w/ Jung ----
# https://github.com/satijalab/seurat/issues/1395


# features to draw heatmap

## RNA assay data  normalization
seurat_integrated.selected.01.combined.name.02.subset@active.assay   # RNA
seurat_integrated.selected.01.combined.name.02.subset <- NormalizeData(seurat_integrated.selected.01.combined.name.02.subset)
head(seurat_integrated.selected.01.combined.name.02.subset@assays$RNA@data)

## Scale RNA assay data
seurat_integrated.selected.01.combined.name.02.subset <- ScaleData(seurat_integrated.selected.01.combined.name.02.subset)
rna.scaled.data_matrix <- seurat_integrated.selected.01.combined.name.02.subset@assays$RNA@scale.data


# gene names that'll be included in the heatmap
# among above gene names, gene names that are overlapping with seurat_integrated.selected.01.combined.name.02.subset.Allmarkers.top5$gene
features_to_draw <- seurat_integrated.selected.01.combined.name.02.subset.Allmarkers.top5$gene
dup.ids <- which(duplicated(features_to_draw))
features_to_draw_no.dup <- features_to_draw[-dup.ids]


## gene names that I'd like to annotate in the heatmap
gene.names <- c("F3",
                "TTYH1",      
                "TIMP2",
                "PTMA",
                "CHD7")

# ## identify the indices of the gene names 
# ids <- which(features_to_draw %in% gene.names)
# 
# ## what to annotate
# genes.to.label <- sample(x = features_to_draw[ids]) 
# labels <- ifelse(features_to_draw %in% gene.names, "black", "transparent")

# labels <- rep(x = "transparent", times = length(x = feature.names.scale.data))
# labels[match(x = genes.to.label, table = feature.names.scale.data)] <- "black"


# complexheatmap
# https://github.com/satijalab/seurat/issues/6214
library(devtools)
install_github("jokergoo/ComplexHeatmap")

browseVignettes("ComplexHeatmap")
## generate a plot
plot_DoHeatmap <- DoHeatmap(seurat_integrated.selected.01.combined.name.02.subset, features = features_to_draw_no.dup,
          slot = "scale.data",
          assay = "RNA",
          disp.max = 4,
          draw.lines = TRUE,lines.width = 40,angle = 0) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") +
  theme(axis.text.y = element_text(color = rev(ifelse(features_to_draw %in% gene.names, "black", "transparent"))))

## >> due to moving/changing indices, the annotation's not done as expected. Other random genes are annotated.
## using ggplot2 geom_text?


plot_DoHeatmap_01 <- DoHeatmap(seurat_integrated.selected.01.combined.name.02.subset, features = features_to_draw_no.dup,
                            slot = "scale.data",
                            assay = "RNA",
                            disp.max = 4,
                            draw.lines = TRUE,lines.width = 40,angle = 0) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")
  
dim(plot_DoHeatmap_01$data)
unique(plot_DoHeatmap_01$data$Feature)

labels <- rep(x = "transparent", times = length(x = nrow(plot_DoHeatmap_01$data)))
labels[which(plot_DoHeatmap_01$data$Feature %in% gene.names)] <- "black"




ifelse(features_to_draw %in% gene.names, "black", "transparent")
unique(plot_DoHeatmap$data$Feature) %in% features_to_draw
features_to_draw %in% unique(plot_DoHeatmap$data$Feature)
which(duplicated(features_to_draw))

plot_DoHeatmap$data[c(1:6, 199:204), ]

plot_DoHeatmap02 <- DoHeatmap(seurat_integrated.selected.01.combined.name.02.subset, features = gene.names,
                            slot = "scale.data",
                            assay = "RNA",
                            disp.max = 4,
                            draw.lines = TRUE,lines.width = 40,angle = 0) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")

all(plot_DoHeatmap$data$Feature %in% features_to_draw)
all(features_to_draw %in% plot_DoHeatmap$data$Feature)
length(features_to_draw)
length(plot_DoHeatmap$data$Feature)

features_to_draw
?DoHeatmap

# + 
#   theme(axis.text.y = element_text(size = 5, color = rev(labels)))

DoHeatmap(seurat_integrated.selected.01.combined.name.02.subset, features = seurat_integrated.selected.01.combined.name.02.subset.Allmarkers.top5$gene,
          slot = "scale.data",
          assay = "SCT",
          disp.max = 4,
          draw.lines = TRUE,lines.width = 40,angle = 0) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + 
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.y = element_text(color = rev(labels)))

rev(labels)


seurat_integrated.selected.01.combined.name.02.subset.Allmarkers.top5$gene[ids]
seurat_integrated.3D <- subset(seurat_integrated.selected.01.combined.name, idents = c("Neuron", "Astrocyte"))
help(VlnPlot)

# Violin Plot ----
seurat_integrated.3D.ast.a5 <- ScaleData(seurat_integrated.3D.ast.a5, verbose = FALSE)
seurat_integrated.3D.ast.a5$tx <- factor(seurat_integrated.3D.ast.a5$tx, levels = c("None",  "Irisin", "Irisin + Anti.Integrin"))
plots <- VlnPlot(seurat_integrated.3D.ast.a5, features = c( "SDCBP", "DDX17"), split.by = "tx", slot = "data", assay = "RNA", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
dev.off()

seurat_integrated.3D.neu.a5 <- ScaleData(seurat_integrated.3D.neu.a5, verbose = FALSE)
seurat_integrated.3D.neu.a5$tx <- factor(seurat_integrated.3D.neu.a5$tx, levels = c("None",  "Irisin", "Irisin + Anti.Integrin"))
plots <- VlnPlot(seurat_integrated.3D.neu.a5, features = c( "HES6", "STMN1"), split.by = "tx", slot = "data", assay = "RNA", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)


View(seurat_integrated.3D@meta.data)

DefaultAssay(seurat_integrated.3D) <- "integrated"
Idents(seurat_integrated.3D) <- "integrated_snn_res.0.4"

DimPlot(seurat_integrated.3D,  label = TRUE, pt.size = 3, label.size = 7) + NoLegend()
seurat_integrated.3D$stim <- factor(seurat_integrated.3D$stim, levels = c("h10_01_b",  "h10_02_b", "h10_03_b"))
DimPlot(seurat_integrated.3D,  label = TRUE, pt.size = 0.8, label.size = 7, split.by = "sample") + NoLegend()
DimPlot(seurat_integrated.3D,  label = TRUE, pt.size = 0.8, label.size = 7, split.by = "day") + NoLegend()
DimPlot(seurat_integrated.3D, label = TRUE, pt.size = 0.8, label.size = 7, group.by = "stim") + NoLegend()

DefaultAssay(seurat_integrated.3D) <- "SCT"
Trimarkers <- c( "GFAP","AQP4", "ID4", "WSB1", "KCNQ1OT1", "SYT1",  "KCNJ6" )
dot01(seurat_integrated.3D, gene = Trimarkers, "Percent expression                 Average gene expression")
dot(seurat_integrated.3D , gene = Trimarkers, "DotPlot (2.0)_by size")

DefaultAssay(seurat_integrated.3D) <- "RNA"



seurat_integrated.3D.Allmarkers <- FindAllMarkers(seurat_integrated.3D, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
seurat_integrated.3D.Allmarkers.top5 <- seurat_integrated.3D.Allmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 40, order_by = avg_log2FC)

View(seurat_integrated.3D.Allmarkers.top5)
View(seurat_integrated.3D.Allmarkers)

DoHeatmap(seurat_integrated.3D, features = seurat_integrated.3D.Allmarkers.top5$gene,
          slot = "scale.data",
          assay = "SCT",
          disp.max = 4,
          draw.lines = TRUE,lines.width = 40,angle = 0
) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 10))


de <- function(x, id1, id2, assay){
  x <- FindMarkers(x,
                   ident.1 = id1,
                   ident.2 = id2,
                   test.use = "wilcox",
                   assay = assay,
                   slot = "data",
                   min.pct = 0.1,
                   logfc.threshold = 0.01)
  return(x)
}


# DE A5. ----
seurat_integrated.3D.ast <- subset(seurat_integrated.3D, idents = c("Astrocyte"))
DefaultAssay(seurat_integrated.3D.ast) <- "RNA"


seurat_integrated.3D.ast.a5 <- subset(seurat_integrated.3D.ast, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b", "a5_03_a", "a5_03_b"))
DimPlot(seurat_integrated.3D.ast.a5)
table(seurat_integrated.3D.ast.a5@meta.data[["sample"]])


# Astrocyte DE by a5_drug ----
Idents(seurat_integrated.3D.ast.a5) <- "tx"
table(seurat_integrated.3D.ast.a5@meta.data[["tx"]])

DefaultAssay(seurat_integrated.3D.ast.a5) <- "RNA"

de.line.list_rna.ast.a5 <- list()
#1.Irisin.vs.None
de.line.list_rna.ast.a5[["Irisin.vs.None"]] <- de(seurat_integrated.3D.ast.a5, 
                                             id1 = "Irisin", 
                                             id2 = "None", 
                                             assay = "RNA")
write.xlsx(de.line.list_rna.ast.a5[["Irisin.vs.None"]], 
           file = "results/seurat_integrated.3D.ast.a5.Irisin.vs.None.xlsx", 
           sheetName = names(de.line.list_rna.ast.a5["Irisin.vs.None"]),
           append = T)



# Count positive/negative logFC within p<=0.2 w/ YJ ----
counting <- function(x){
  library(dplyr)
  p.filter <- x %>% filter(p_val_adj <= 0.2)
  positive.FC <- nrow(p.filter %>% filter(avg_log2FC > 0))
  negative.FC <- nrow(p.filter %>% filter(avg_log2FC < 0))
  print(paste0("pos = ", positive.FC, " ", "neg = ", negative.FC))
}

lapply(de.line.list_rna.ast.a5, counting)

# :: Ast.GO.Irisin.vs.None ----
# rm(Astrocyte.GO.Irisin.vs.None)
Astrocyte.GO.Irisin.vs.None <- de.line.list_rna.ast.a5[["Irisin.vs.None"]] 
Astrocyte.GO.Irisin.vs.None$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and p_val < 0.05, set as "UP" 
Astrocyte.GO.Irisin.vs.None.up <- Astrocyte.GO.Irisin.vs.None$diffexpressed[Astrocyte.GO.Irisin.vs.None$avg_log2FC > 0.2 & Astrocyte.GO.Irisin.vs.None$p_val_adj < 0.2] <- "UP"
Astrocyte.GO.Irisin.vs.None.up.ids <- which( Astrocyte.GO.Irisin.vs.None$diffexpressed == "UP") 
Astrocyte.GO.Irisin.vs.None.up <- Astrocyte.GO.Irisin.vs.None[Astrocyte.GO.Irisin.vs.None.up.ids, ]
Astrocyte.GO.Irisin.vs.None.genes_to_test.up <- rownames(Astrocyte.GO.Irisin.vs.None.up)

help(enrichGO)
# genes_to_test <- "TREM2"
Astrocyte.GO.Irisin.vs.None.up.GO_results <- enrichGO(gene = Astrocyte.GO.Irisin.vs.None.genes_to_test.up, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont ="BP")
# help("enrichGO")

# fit <- plot(barplot(GO_results, showCategory = 20))
plot(dotplot(Astrocyte.GO.Irisin.vs.None.up.GO_results, showCategory = 15))
cnetplot(Astrocyte.GO.Irisin.vs.None.up.GO_results)



#2.Irisin.vs.Irisin + Anti.Integrin
de.line.list_rna.ast.a5[["Irisin.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.ast.a5, 
                                                      id1 = "Irisin", 
                                                      id2 = "Irisin + Anti.Integrin", 
                                                      assay = "RNA")
write.xlsx(de.line.list_rna.ast.a5[["Irisin.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.ast.a5.Irisin.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.ast.a5["Irisin.vs.Irisin + Anti.Integrin"]),
           append = T)

#3.None.vs.Irisin + Anti.Integrin
de.line.list_rna.ast.a5[["None.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.ast.a5, 
                                                        id1 = "None", 
                                                        id2 = "Irisin + Anti.Integrin", 
                                                        assay = "RNA")
write.xlsx(de.line.list_rna.ast.a5[["None.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.ast.a5.None.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.ast.a5["None.vs.Irisin + Anti.Integrin"]),
           append = T)

# # export as xlsx
# library(xlsx)
# for(i in 1:length(de.line.list_rna.ast)){
#   write.xlsx(de.line.list_rna.ast[[i]], 
#              file = "results/seurat_integrated.3D.ast.a5.02.vs.01.xlsx", 
#              sheetName = names(de.line.list_rna.ast[i]),
#              append = T)
# }

seurat_integrated.3D.neu <- subset(seurat_integrated.3D, idents = c("Neuron"))
table(seurat_integrated.3D.neu@meta.data[["sample"]])
Idents(seurat_integrated.3D.neu) <- "sample"
seurat_integrated.3D.neu.a5 <- subset(seurat_integrated.3D.neu, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b", "a5_03_a", "a5_03_b"))

# Neuron DE by a5_drug ----
table(seurat_integrated.3D.neu.a5@meta.data[["tx"]])
Idents(seurat_integrated.3D.neu.a5) <- "tx"

de.line.list_rna.neu.a5 <- list()
#1.Tx1.vs.Ctrl
de.line.list_rna.neu.a5[["Irisin.vs.None"]] <- de(seurat_integrated.3D.neu.a5, 
                                               id1 = "Irisin", 
                                               id2 = "None", 
                                               assay = "RNA")
write.xlsx(de.line.list_rna.neu.a5[["Irisin.vs.None"]], 
           file = "results/seurat_integrated.3D.neu.a5.Irisin.vs.None.xlsx", 
           sheetName = names(de.line.list_rna.neu.a5["Irisin.vs.None"]),
           append = T)


# :: Neuron.GO.Irisin.vs.None ----
Neuron.GO.Irisin.vs.None <- de.line.list_rna.neu.a5[["Irisin.vs.None"]] 
Neuron.GO.Irisin.vs.None$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and p_val < 0.05, set as "UP" 
Neuron.GO.Irisin.vs.None.up <- Neuron.GO.Irisin.vs.None$diffexpressed[Neuron.GO.Irisin.vs.None$avg_log2FC > 0.1 & Neuron.GO.Irisin.vs.None$p_val_adj < 0.2] <- "UP"
Neuron.GO.Irisin.vs.None.up.ids <- which( Neuron.GO.Irisin.vs.None$diffexpressed == "UP") 
Neuron.GO.Irisin.vs.None.up <- Neuron.GO.Irisin.vs.None[Neuron.GO.Irisin.vs.None.up.ids, ]
Neuron.GO.Irisin.vs.None.up.genes_to_test.up <- rownames(Neuron.GO.Irisin.vs.None.up)

# genes_to_test <- "TREM2"
Neuron.GO.Irisin.vs.None.up.GO_results <- enrichGO(gene = Neuron.GO.Irisin.vs.None.up.genes_to_test.up, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont ="ALL")
# help("enrichGO")

# fit <- plot(barplot(GO_results, showCategory = 20))
plot(dotplot(Neuron.GO.Irisin.vs.None.up.GO_results, showCategory = 15))
cnetplot(Neuron.GO.Irisin.vs.None.up.GO_results)


#2.Irisin.vs.Irisin + Anti.Integrin
de.line.list_rna.neu.a5[["Irisin.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.neu.h10, 
                                                        id1 = "Irisin", 
                                                        id2 = "Irisin + Anti.Integrin", 
                                                        assay = "RNA")
write.xlsx(de.line.list_rna.neu.a5[["Irisin.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.neu.a5.Irisin.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.neu.a5["Irisin.vs.Irisin + Anti.Integrin"]),
           append = T)
#3.None.vs.Irisin + Anti.Integrin
de.line.list_rna.neu.a5[["None.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.neu.h10, 
                                                         id1 = "None", 
                                                         id2 = "Irisin + Anti.Integrin", 
                                                         assay = "RNA")
write.xlsx(de.line.list_rna.neu.a5[["None.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.neu.a5.None.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.neu.a5["None.vs.Irisin + Anti.Integrin"]),
           append = T)

FeaturePlot(seurat_integrated.3D, pt.size = 2.0, features = c("GFAP", "AQP4"))
FeaturePlot(seurat_integrated.3D, pt.size = 2.0, features = c("WSB1", "SYT1", "KCNQ1OT1"))


# DE H10. ----
seurat_integrated.3D.ast.h10@meta.data[["tx"]] <- factor(seurat_integrated.3D.ast.h10@meta.data[["tx"]], levels = c( "None", "Irisin", "Irisin + Anti.Integrin"))

seurat_integrated.3D.ast.h10@m
VlnPlot(seurat_integrated.3D.ast.h10, features = c( "GFAP", "VIM"), split.by = "tx", slot = "data", assay = "RNA", y.max = 100,
                 pt.size = 0, combine = FALSE)
plots <- VlnPlot(seurat_integrated.3D.ast.h10, features = c( "GFAP", "VIM"), split.by = "tx", slot = "data", assay = "RNA", y.max = 50,
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

Idents(seurat_integrated.3D) <- "sample"
seurat_integrated.3D.h10 <- subset(seurat_integrated.3D, idents = c("h10_01_a", "h10_01_b", "h10_02_b", "h10_03_a", "h10_03_b" ))
Idents(seurat_integrated.3D.h10) <- "SCT_snn_res.0.4"
seurat_integrated.3D.h10.name <- RenameIdents(object = seurat_integrated.3D.h10, 
                                                            "0" = "Astrocyte",
                                                            "1" = "Astrocyte",
                                                            
                                                            "3" =  "Neuron",
                                                            "4" = "Neuron",
                                                            "5" =  "Astrocyte",
                                                            
                                                            "7" = "Astrocyte")


DimPlot(seurat_integrated.3D.h10.name, label = TRUE, pt.size = 3, label.size = 7) + NoLegend()
DefaultAssay(seurat_integrated.3D.h10.name) <- "SCT"
Trimarkers <- c( "GFAP",   "AQP4", "ID4",   "WSB1", "KCNQ1OT1", "SYT1", "KCNJ6" )
dot01(seurat_integrated.3D.h10.name, gene = Trimarkers, "Percent expression      Average gene expression")
dot(seurat_integrated.3D.h10.name , gene = Trimarkers, "DotPlot (2.0)_by size")



table(seurat_integrated.3D.ast@meta.data[["sample"]])
Idents(seurat_integrated.3D.ast) <- "sample"
seurat_integrated.3D.ast.h10 <- subset(seurat_integrated.3D.ast, idents = c("h10_01_a", "h10_01_b", "h10_02_b", "h10_03_a", "h10_03_b" ))


# Astrocyte DE by h10_drug ----
Idents(seurat_integrated.3D.ast.h10) <- "tx"
table(seurat_integrated.3D.ast.h10@meta.data[["tx"]])
dev.off()
DefaultAssay(seurat_integrated.3D.ast.h10) <- "RNA"

de.line.list_rna.ast.h10 <- list()
#1.Irisin.vs.None
de.line.list_rna.ast.h10[["Irisin.vs.None"]] <- de(seurat_integrated.3D.ast.h10, 
                                                  id1 = "Irisin", 
                                                  id2 = "None", 
                                                  assay = "RNA")
library(xlsx)
write.xlsx(de.line.list_rna.ast.h10[["Irisin.vs.None"]], 
           file = "results/seurat_integrated.3D.ast.h10.Irisin.vs.None.xlsx", 
           sheetName = names(de.line.list_rna.ast.h10["Irisin.vs.None"]),
           append = T)

#2.Irisin.vs.Irisin + Anti.Integrin
de.line.list_rna.ast.h10[["Irisin.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.ast.h10, 
                                                                    id1 = "Irisin", 
                                                                    id2 = "Irisin + Anti.Integrin", 
                                                                    assay = "RNA")
write.xlsx(de.line.list_rna.ast.h10[["Irisin.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.ast.h10.Irisin.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.ast.h10["Irisin.vs.Irisin + Anti.Integrin"]),
           append = T)

#3.None.vs.Irisin + Anti.Integrin
de.line.list_rna.ast.h10[["None.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.ast.h10, 
                                                                  id1 = "None", 
                                                                  id2 = "Irisin + Anti.Integrin", 
                                                                  assay = "RNA")
write.xlsx(de.line.list_rna.ast.h10[["None.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.ast.h10.None.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.ast.h10["None.vs.Irisin + Anti.Integrin"]),
           append = T)

# # export as xlsx
# library(xlsx)
# for(i in 1:length(de.line.list_rna.ast)){
#   write.xlsx(de.line.list_rna.ast[[i]], 
#              file = "results/seurat_integrated.3D.ast.a5.02.vs.01.xlsx", 
#              sheetName = names(de.line.list_rna.ast[i]),
#              append = T)
# }
# Neuron DE by h10_drug ----
seurat_integrated.3D.neu <- subset(seurat_integrated.3D, idents = c("Neuron"))
table(seurat_integrated.3D.neu@meta.data[["sample"]])
Idents(seurat_integrated.3D.neu) <- "sample"
seurat_integrated.3D.neu.h10 <- subset(seurat_integrated.3D.neu, idents = c("h10_01_a", "h10_01_b", "h10_02_b", "h10_03_a", "h10_03_b"))


table(seurat_integrated.3D.neu.h10@meta.data[["tx"]])
Idents(seurat_integrated.3D.neu.h10) <- "tx"

de.line.list_rna.neu.h10 <- list()
#1.Tx1.vs.Ctrl
de.line.list_rna.neu.h10[["Irisin.vs.None"]] <- de(seurat_integrated.3D.neu.h10, 
                                                  id1 = "Irisin", 
                                                  id2 = "None", 
                                                  assay = "RNA")
write.xlsx(de.line.list_rna.neu.h10[["Irisin.vs.None"]], 
           file = "results/seurat_integrated.3D.neu.h10.Irisin.vs.None.xlsx", 
           sheetName = names(de.line.list_rna.neu.h10["Irisin.vs.None"]),
           append = T)

#2.Irisin.vs.Irisin + Anti.Integrin
de.line.list_rna.neu.h10[["Irisin.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.neu.h10, 
                                                                    id1 = "Irisin", 
                                                                    id2 = "Irisin + Anti.Integrin", 
                                                                    assay = "RNA")
write.xlsx(de.line.list_rna.neu.h10[["Irisin.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.neu.h10.Irisin.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.neu.h10["Irisin.vs.Irisin + Anti.Integrin"]),
           append = T)
#3.None.vs.Irisin + Anti.Integrin
de.line.list_rna.neu.h10[["None.vs.Irisin + Anti.Integrin"]] <- de(seurat_integrated.3D.neu.h10, 
                                                                  id1 = "None", 
                                                                  id2 = "Irisin + Anti.Integrin", 
                                                                  assay = "RNA")
write.xlsx(de.line.list_rna.neu.h10[["None.vs.Irisin + Anti.Integrin"]], 
           file = "results/seurat_integrated.3D.neu.h10.None.vs.Irisin + Anti.Integrin.xlsx", 
           sheetName = names(de.line.list_rna.neu.h10["None.vs.Irisin + Anti.Integrin"]),
           append = T)


# Heatmap.A5.Ast.Irisin.vs.None ----
seurat_integrated.3D.ast <- subset(seurat_integrated.3D, idents = c("Astrocyte"))

table(seurat_integrated.3D.ast.a5@meta.data[["sample"]])
Idents(seurat_integrated.3D.ast) <- "sample"
seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None <- subset(seurat_integrated.3D.ast, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b"))

DoHeatmap(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, 
          features = c("STEAP4","S1PR3","TIMP1","HSPB1","CD44","OSMR","CP","VIM","GFAP",                                                                                  
                       "SERPING","FBLN5","UGT1A1","FKBP5","PSMB8","AMIGO2",
                       "PTX3","S100A10","CD109","EMP1","B3GNT5"),
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 10))


DoHeatmap(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, features = c("GFAP",
                                                                                  "ALDOC",
                                                                                  "C3",
                                                                                  "FABP7",
                                                                                  "CHI3L1",
                                                                                  "CRYAB",
                                                                                  "BLBP",
                                                                                  "HSPB1",
                                                                                  "HSP27",                                                                                  
                                                                                  "IL17R",
                                                                                  "LCN2",
                                                                                  "MAO-B",
                                                                                  "NTRK2",
                                                                                  "TRKB",
                                                                                  "S100B",
                                                                                  "SERPINA3N",
                                                                                  "EAAT2",
                                                                                  "EAAT1",
                                                                                  "TSPO",
                                                                                  "SOX9",
                                                                                  "STAT3",
                                                                                  "SYNEMIN",
                                                                                  "THBS-1",
                                                                                  "ACT",
                                                                                  "NESTIN",
                                                                                  "NFAT",
                                                                                  "VIMENTIN",
                                                                                  "YKL40",
                                                                                  "KIR4.1",
                                                                                  
                                                                                  "VIM",
                                                                                  
                                                                                  "ASPG",
                                                                                  "CH25H",
                                                                                  "ZBP1",
                                                                                  "TLR2",
                                                                                  "A2M",
                                                                                  "CXCL16",
                                                                                  "F830016B08RIK",
                                                                                  "GBP4",
                                                                                  "H2-D1",
                                                                                  "H2-Q4",
                                                                                  "H2-Q6",
                                                                                  "HCAR2",
                                                                                  "ICAM1",
                                                                                  "IER3",
                                                                                  "IFI47",
                                                                                  "IFIT1",
                                                                                  "IIGP1",
                                                                                  "IL1A",
                                                                                  "IL1RN",
                                                                                  "NLRC5",
                                                                                  "TAP1",
                                                                                  "TRIM30A","CXCL10","DTX3L","FKBP5","GBP6","IFITM3","CP"),
          
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))

seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None <- subset(seurat_integrated.3D.ast, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b", "a5_03_a",  "a5_03_b"))
Idents(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None) <- "tx"
DefaultAssay(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None) <- "RNA"
seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None <- NormalizeData(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, verbose = TRUE)
seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None <- ScaleData(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, verbose = FALSE)

seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01 <- AverageExpression(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, return.seurat = TRUE, assays = "RNA", slot = "scale.data")
write.xlsx(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@assays[["RNA"]]@scale.data,
           file = "results/3D.ast.a5.heatmap.Irisin.vs.None.ave.01.xlsx")

seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@active.ident <- factor(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@active.ident, levels = c("Irisin", "None"))
DoHeatmap(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01, features = c("GFAP",
                                                                                  "ALDOC",
                                                                                  "C3"),
          
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))

# :: Heatmap.A5.Ast.Gene.name.label ----

genes.to.label <- c("PTN", "GFAP")
genes.to.label <- sample(x = seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@assays[["RNA"]]@scale.data$gene[1:100]) 
labels <- rep(x = "transparent", times = length(x = seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@assays[["RNA"]]@scale.data$gene))
labels[match(x = genes.to.label, table = seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01@assays[["RNA"]]@scale.data$gene)] <- "black"

help("DoHeatmap")
DoHeatmap(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None.ave.01, features = c("GFAP",
                                                                                  "ALDOC",
                                                                                  "C3"),
          
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))


# Astrocyte DE by a5_drug ----



# ave.Heatmap.A5.Neu.Irisin.vs.None ----
seurat_integrated.3D.neu <- subset(seurat_integrated.3D, idents = c("Neuron"))

table(seurat_integrated.3D.neu@meta.data[["sample"]])
Idents(seurat_integrated.3D.neu) <- "sample"
seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None <- subset(seurat_integrated.3D.neu, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b", "a5_03_a",  "a5_03_b"))


Idents(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None) <- "tx"
DefaultAssay(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None) <- "RNA"
seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None <- NormalizeData(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None, verbose = TRUE)
seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None <- ScaleData(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None, verbose = FALSE)

seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None.ave.01 <- AverageExpression(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None, return.seurat = TRUE, assays = "RNA", slot = "scale.data")
write.xlsx(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None.ave.01@assays[["RNA"]]@scale.data,
           file = "results/3D.neu.a5.heatmap.ave.01.xlsx")

seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None.ave.01@active.ident <- factor(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None.ave.01@active.ident, levels = c("Irisin", "None"))
DoHeatmap(seurat_integrated.3D.neu.a5.heatmap.Irisin.vs.None.ave.01, features = c('GADD45G', "HNRNPA1",
                                                                                  "GLCCI1", "GAS5", "BTG1", "AFDN",
                                                                                  "KCNQ1OT1", "STMN1",  'KIF21A','CD44', 'CDKN1A', 'CST3', 'ATP1B2', 'CLU', 'GABBR2', 'DCLK1', 'CADM1', 'CDKN1A', 'APOE', "PLPP3", "CD44", "SLC6A11", "GPX3"),
          
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))

# ave.Heatmap.h10.Ast.tx ----
seurat_integrated.3D.ast <- subset(seurat_integrated.3D, idents = c("Astrocyte"))

features <- c("RELA", "GFAP" , "BDNF", "MME", "NFKB1", "NFKB2")
VlnPlot(seurat_integrated.3D.ast.a5.heatmap.Irisin.vs.None, features = features, slot = "counts", log = TRUE)

table(seurat_integrated.3D.ast.h10@meta.data[["sample"]])
Idents(seurat_integrated.3D.ast.h10) <- "sample"
seurat_integrated.3D.ast.h10.heatmap <- subset(seurat_integrated.3D.ast.h10, idents = c("h10_01_a", "h10_01_b", "h10_02_b"))
Idents(seurat_integrated.3D.ast.h10.heatmap) <- "tx"
table(seurat_integrated.3D.ast.h10.heatmap@meta.data[["tx"]])
DefaultAssay(seurat_integrated.3D.ast.h10.heatmap) <- "RNA"
seurat_integrated.3D.ast.h10.heatmap <- NormalizeData(seurat_integrated.3D.ast.h10.heatmap, verbose = TRUE)
seurat_integrated.3D.ast.h10.heatmap <- ScaleData(seurat_integrated.3D.ast.h10.heatmap, verbose = FALSE)

seurat_integrated.3D.ast.h10.heatmap.ave.01 <- AverageExpression(seurat_integrated.3D.ast.h10.heatmap, return.seurat = TRUE)
write.xlsx(seurat_integrated.3D.ast.h10.heatmap.ave.01@assays[["RNA"]]@scale.data,
           file = "results/3D.ast.h10.heatmap.ave.01_Nov.23.xlsx")

seurat_integrated.3D.ast.h10.heatmap.ave.01@meta.data[["tx"]] <- factor(seurat_integrated.3D.ast.h10.heatmap.ave.01@meta.data[["tx"]], levels = c( "None", "Irisin", "Irisin + Anti.Integrin"))

DoHeatmap(seurat_integrated.3D.ast.h10.heatmap.ave.01, 
          features = c("STEAP4","S1PR3","TIMP1","HSPB1","CD44","OSMR","CP","VIM","GFAP",                                                                                  
                       "SERPING","FBLN5","UGT1A1","FKBP5","PSMB8","AMIGO2",
                       "PTX3","S100A10","CD109","EMP1","B3GNT5"),
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = NULL,
          disp.max = NULL,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))

# : Reactive heatmap overlap ----
reactive <- read.table("results/reactive.astrocyte.txt")
row.names(reactive ) <- reactive$V1

# DAM.down.ChR2 <- read.table("results/DAM_DOWN_1.ChR2 vs GFP.txt", header=TRUE)
# row.names(DAM.down.ChR2) <- DAM.down.ChR2$Row.names
# 
# DAM.up.Genus <- read.table("results/DAM_UP_2.Genus vs GFP.txt", header=TRUE)
# row.names(DAM.up.Genus) <- DAM.up.Genus$Row.names
# 
# DAM.down.Genus <- read.table("results/DAM_DOWN_2.Genus vs GFP.txt", header=TRUE)
# row.names(DAM.down.Genus) <- DAM.down.Genus$Row.names


Reactive_overlap_values <- list()
Reactive_overlap_values[[1]] <- merge(reactive, seurat_integrated.3D.ast.h10.heatmap.ave.01@assays[["RNA"]]@scale.data,  by = "row.names")
Reactive_overlap_values[[2]] <- merge(reactive, seurat_integrated.3D.ast.h10.heatmap.ave.01@assays[["RNA"]]@data,  by = "row.names")

# : ave.Heatmap. h10 separated ----
Reactive_overlap_values[[1]]
row.names(Reactive_overlap_values[[1]]) <- Reactive_overlap_values[[1]]$Row.names
Reactive_overlap_values[[1]]$Row.names <-NULL
Reactive_overlap_values[[1]]$V1 <-NULL

write.xlsx(Reactive_overlap_values[[1]],
           file = "results/3D.ast.h10.heatmap.ave.reactive.xlsx")
reactive <- read.table("results/reactive.astrocyte.txt")


ave.heat <- read.table("results/reactive.astrocyte.ordered.txt")
ave.heat.data.mat <- data.matrix(ave.heat)
ave.heat.data.mat.scale <- t(scale(t(ave.heat.data.mat)))
# seurat_integrated.3D.neu.h10.heatmap.ave.sep <- as.data.frame(AverageExpression(seurat_integrated.3D.neu.h10.heatmap,  assays = "RNA", slot = "data"))
# seurat_integrated.3D.neu.h10.heatmap.ave.sep <- data.matrix(seurat_integrated.3D.neu.h10.heatmap.ave.sep)

Heatmap(ave.heat.data.mat.scale, name = "Reactive Astrocyte Genes in H10", 
        row_order = NULL,
        column_order = NULL)

Heatmap(ave.heat.data.mat.scale,  name = "Normalized Heatmap: Reactive astrocyte genes",
        # row_order = NULL,
        # row_order = sort(rownames(ave.heat.data.mat.scale)), 
        column_order = sort(colnames(ave.heat.data.mat.scale)),
        # cluster_row_slices = FALSE, 
        # cluster_column_slices = FALSE,
        column_km = 1,
        
        row_gap = unit(2, "mm"),
        column_gap = unit(c(7, 7), "mm"),
        # border = TRUE,
        row_split = factor(c(rep("PAN", times = 10), rep("A1", times = 7), rep("A2", times = 7)), levels = c("PAN", "A1", "A2")),
        column_split = factor(c(rep("None", times = 1), rep("Irisin", times = 1), rep("Irisin + Anti.Integrin", times = 1)), levels = c("None", "Irisin", "Irisin + Anti.Integrin")))
# cluster_rows = FALSE,)

# ave.Heatmap.h10.neu.tx ----
seurat_integrated.3D.ast <- subset(seurat_integrated.3D, idents = c("Astrocyte"))

table(seurat_integrated.3D.neu.h10@meta.data[["sample"]])
Idents(seurat_integrated.3D.neu.h10) <- "sample"
seurat_integrated.3D.neu.h10.heatmap <- subset(seurat_integrated.3D.neu.h10, idents = c("h10_01_a", "h10_01_b", "h10_02_b", "h10_03_a", "h10_03_b"))
Idents(seurat_integrated.3D.neu.h10.heatmap) <- "tx"
DefaultAssay(seurat_integrated.3D.neu.h10.heatmap) <- "RNA"
seurat_integrated.3D.neu.h10.heatmap <- NormalizeData(seurat_integrated.3D.neu.h10.heatmap, verbose = TRUE)
seurat_integrated.3D.neu.h10.heatmap <- ScaleData(seurat_integrated.3D.neu.h10.heatmap, verbose = FALSE)

seurat_integrated.3D.neu.h10.heatmap.ave.01 <- AverageExpression(seurat_integrated.3D.neu.h10.heatmap, return.seurat = TRUE, assays = "RNA", slot = "scale.data")
write.xlsx(seurat_integrated.3D.neu.h10.heatmap.ave.01@assays[["RNA"]]@scale.data,
           file = "results/seurat_integrated.3D.neu.h10.heatmap.ave.01.xlsx")

seurat_integrated.3D.neu.h10.heatmap.ave.01@meta.data[["tx"]] <- factor(seurat_integrated.3D.neu.h10.heatmap.ave.01@meta.data[["tx"]], levels = c( "None", "Irisin", "Irisin + Anti.Integrin"))



DoHeatmap(seurat_integrated.3D.neu.h10.heatmap.ave.01, features = c("STEAP4","S1PR3","TIMP1","HSPB1","CD44","OSMR","CP","VIM","GFAP",                                                                                  
                                                                    "SERPING","FBLN5","UGT1A1","FKBP5","PSMB8","AMIGO2",
                                                                    "PTX3","S100A10","CD109","EMP1","B3GNT5"),
          
          group.bar = TRUE,
          group.colors = NULL,
          disp.min = -0.2,
          disp.max = 0.2,
          slot = "scale.data",
          assay = "RNA",
          draw.lines = TRUE, 
          lines.width = NULL,
          angle = 0,
          label = TRUE,
          size = 5,
          hjust = 0.5,
          group.bar.height = 0.01,
          combine = TRUE
)  + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.y = element_text(size = 15))

# ave.Heatmap. h10 separated ----
ave.heat <- read.table("results/reactive.astrocyte.ordered.txt", head = TRUE)
ave.heat.data.mat <- data.matrix(ave.heat)
ave.heat.data.mat.scale <- t(scale(t(ave.heat.data.mat)))
# seurat_integrated.3D.neu.h10.heatmap.ave.sep <- as.data.frame(AverageExpression(seurat_integrated.3D.neu.h10.heatmap,  assays = "RNA", slot = "data"))
# seurat_integrated.3D.neu.h10.heatmap.ave.sep <- data.matrix(seurat_integrated.3D.neu.h10.heatmap.ave.sep)

Heatmap(ave.heat.data.mat.scale, name = "Reactive Astrocyte Genes in H10", row_order = NULL,
        column_order = NULL)

Heatmap(ave.heat.data.mat.scale,  name = "Normalized Heatmap: Reactive astrocyte genes",
        row_order = NULL,
        # row_order = sort(rownames(ave.heat.data.mat.scale)), 
        column_order = sort(colnames(ave.heat.data.mat.scale)),
        row_gap = unit(2, "mm"),
        column_gap = unit(c(7, 7), "mm"),
        # border = TRUE,
        row_split = factor(c(rep("PAN", times = 10), rep("A1", times = 7), rep("A2", times = 7)), levels = c("PAN", "A1", "A2")),
        column_split = factor(c(rep("None", times = 1), rep("Irisin", times = 1), rep("Irisin + Anti.Integrin", times = 1)), levels = c("None", "Irisin", "Irisin + Anti.Integrin")))
        # cluster_rows = FALSE,)


dim(ave.heat)


scaled_mat.ave.heat = t(scale(t(ave.heat)))




# Export for Vln.Plot ----


table(seurat_integrated.selected.01.combined.name$sample)

Idents(seurat_integrated.3D.ast) <- "sample"
table(Idents(seurat_integrated.selected.01.combined.name))

rm(seurat_Ast_A5.AB)
seurat_Ast_A5.AB <- subset(seurat_integrated.3D.ast, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b",  "a5_03_a",  "a5_03_b"))
DimPlot(seurat_Ast_A5.AB,  label = TRUE, pt.size = 3, label.size = 7) + NoLegend()

seurat_Ast_A5.AB$sample <- factor(seurat_Ast_A5.AB$sample, levels = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b",  "a5_03_a",  "a5_03_b"))

# A5 ast vln export ----
DimPlot(seurat_integrated.3D.ast.a5,   pt.size = 0.8, label.size = 7, split.by = "tx") + NoLegend()
seurat_integrated.3D.ast.a5.None <- subset(seurat_integrated.3D.ast.a5, idents = "None")
write.table(seurat_integrated.3D.ast.a5.None@assays[["RNA"]]@data, file='seurat_integrated.3D.ast.a5.None.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.a5.Irisin <- subset(seurat_integrated.3D.ast.a5, idents = "Irisin")
write.table(seurat_integrated.3D.ast.a5.Irisin@assays[["RNA"]]@data, file='seurat_integrated.3D.ast.a5.Irisin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.a5.Irisin.Anti.Integrin <- subset(seurat_integrated.3D.ast.a5, idents = "Irisin + Anti.Integrin")
write.table(seurat_integrated.3D.ast.a5.Irisin.Anti.Integrin@assays[["RNA"]]@data, file='seurat_integrated.3D.ast.a5.Irisin + Anti.Integrin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# H10 neu vln export ----
Idents(seurat_integrated.3D.neu.h10) <- "tx"
seurat_integrated.3D.neu.h10.None <- subset(seurat_integrated.3D.neu.h10, idents = "None")
write.table(seurat_integrated.3D.neu.h10.None@assays[["RNA"]]@data, file='3D.neu.h10.None.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.neu.h10.Irisin <- subset(seurat_integrated.3D.neu.h10, idents = "Irisin")
write.table(seurat_integrated.3D.neu.h10.Irisin@assays[["RNA"]]@data, file='3D.neu.h10.Irisin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.neu.h10.Irisin.Anti.Integrin <- subset(seurat_integrated.3D.neu.h10, idents = "Irisin + Anti.Integrin")
write.table(seurat_integrated.3D.neu.h10.Irisin.Anti.Integrin@assays[["RNA"]]@data, file='3D.neu.h10.Irisin + Anti.Integrin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# H10 ast vln export ----
Idents(seurat_integrated.3D.ast.h10) <- "tx"
seurat_integrated.3D.ast.h10.None <- subset(seurat_integrated.3D.ast.h10, idents = "None")
write.table(seurat_integrated.3D.ast.h10.None@assays[["RNA"]]@data, file='3D.ast.h10.None.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.h10.Irisin <- subset(seurat_integrated.3D.ast.h10, idents = "Irisin")
write.table(seurat_integrated.3D.ast.h10.Irisin@assays[["RNA"]]@data, file='3D.ast.h10.Irisin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.h10.Irisin.Anti.Integrin <- subset(seurat_integrated.3D.ast.h10, idents = "Irisin + Anti.Integrin")
write.table(seurat_integrated.3D.ast.h10.Irisin.Anti.Integrin@assays[["RNA"]]@data, file='3D.ast.h10.Irisin + Anti.Integrin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# a5 neu vln export ----
Idents(seurat_integrated.3D.neu.a5) <- "tx"
seurat_integrated.3D.neu.a5.None <- subset(seurat_integrated.3D.neu.a5, idents = "None")
write.table(seurat_integrated.3D.neu.a5.None@assays[["RNA"]]@data, file='3D.neu.a5.None.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.neu.a5.Irisin <- subset(seurat_integrated.3D.neu.a5, idents = "Irisin")
write.table(seurat_integrated.3D.neu.a5.Irisin@assays[["RNA"]]@data, file='3D.neu.a5.Irisin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.neu.a5.Irisin.Anti.Integrin <- subset(seurat_integrated.3D.neu.a5, idents = "Irisin + Anti.Integrin")
write.table(seurat_integrated.3D.neu.a5.Irisin.Anti.Integrin@assays[["RNA"]]@data, file='3D.neu.a5.Irisin + Anti.Integrin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

# a5 ast vln export ----
Idents(seurat_integrated.3D.ast.a5) <- "tx"
seurat_integrated.3D.ast.a5.None <- subset(seurat_integrated.3D.ast.a5, idents = "None")
write.table(seurat_integrated.3D.ast.a5.None@assays[["RNA"]]@data, file='3D.ast.a5.None.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.a5.Irisin <- subset(seurat_integrated.3D.ast.a5, idents = "Irisin")
write.table(seurat_integrated.3D.ast.a5.Irisin@assays[["RNA"]]@data, file='3D.ast.a5.Irisin.tsv', quote=FALSE, sep='\t', col.names = TRUE)

seurat_integrated.3D.ast.a5.Irisin.Anti.Integrin <- subset(seurat_integrated.3D.ast.a5, idents = "Irisin + Anti.Integrin")
write.table(seurat_integrated.3D.ast.a5.Irisin.Anti.Integrin@assays[["RNA"]]@data, file='3D.ast.a5.Irisin + Anti.Integrin.tsv', quote=FALSE, sep='\t', col.names = TRUE)




table(seurat_integrated.selected.01.combined@meta.data$integrated_snn_res.0.4, seurat_integrated.selected.01.combined@meta.data$sample)

rm(seurat_integrated.3D.neu.a5.02.vs.01)
table(seurat_integrated.3D.neu@meta.data[["sample"]])
# a5_Neuron
seurat_integrated.3D.neu <- subset(seurat_integrated.3D, idents = c("Neuron"))
Idents(seurat_integrated.3D.neu) <- "sample"
seurat_integrated.3D.neu.a5.02.vs.01 <- subset(seurat_integrated.3D.neu, idents = c("a5_01_a", "a5_01_b", "a5_02_a",  "a5_02_b"))
DimPlot(seurat_integrated.3D.neu.a5.02.vs.01)
# Neuron DE by line
Idents(seurat_integrated.3D.neu.a5.02.vs.01) <- "sample"

de.line.list_rna <- list()
de.line.list_rna[["a5_02_a vs. a5_01_a"]] <- de(seurat_integrated.3D.neu.a5.02.vs.01, 
                                                id1 = "a5_02_a", 
                                                id2 = "a5_01_a", 
                                                assay = "RNA")
de.line.list_rna[["a5_02_b vs. a5_01_b"]] <- de(seurat_integrated.3D.neu.a5.02.vs.01, 
                                                id1 = "a5_02_b", 
                                                id2 = "a5_01_b", 
                                                assay = "RNA")
table(seurat_integrated.3D@meta.data[["stim"]])
# export as xlsx
library(xlsx)
for(i in 1:length(de.line.list_rna)){
  write.xlsx(de.line.list_rna[[i]], 
             file = "results/seurat_integrated.3D.neu.a5.02.vs.01.xlsx", 
             sheetName = names(de.line.list_rna[i]),
             append = T)
}


#&&&&&&& OLD Psudobulk &&&&&&& ----
View(seurat_integrated.3D.ast.a5.02.vs.01@meta.data)
seurat_integrated.3D.ast.a5.02.vs.01$pseudo <- paste0(seurat_integrated.3D.ast.a5.02.vs.01$sample, seurat_integrated.3D.ast.a5.02.vs.01$batch)

DefaultAssay(seurat_integrated.3D.ast.a5.02.vs.01)

cts <- AggregateExpression(seurat_integrated.3D.ast.a5.02.vs.01, 
                           group.by = c("sample", "pseudo"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA
view(cts)
# # transpose
# cts.t <- t(cts)
# 
# 
# # convert to data.frame
# cts.t <- as.data.frame(cts.t)
# 
# # get values where to split
# splitRows <- gsub('_.*', '', rownames(cts.t))
# 
# 
# # split data.frame
# cts.split <- split.data.frame(cts.t,
#                               f = factor(splitRows))
# 
# # fix colnames and transpose
# 
# cts.split.modified <- lapply(cts.split, function(x){
#   rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
#   t(x)
#   
# })
# 
# #gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Let's run DE analysis with B cells
# 1. Get counts matrix

colnames(as.data.frame(cts))

# counts_astro.neuron <- list()
# counts_astro.neuron[["astro_h10_01.vsh10_02"]] <- as.data.frame(cts[, c("Astrocyte_h10_01_bCtrl", "Astrocyte_h10_02_bTx1")])
# counts_astro.neuron[["astro_h10_02.vsh10_03"]] <- as.data.frame(cts[, c("Astrocyte_h10_02_bTx1", "Astrocyte_h10_03_bIGG")])
# counts_astro.neuron[["neuron_h10_01.vsh10_02"]] <- as.data.frame(cts[, c("Neuron_h10_01_bCtrl", "Neuron_h10_02_bTx1")])
# counts_astro.neuron[["neuron_h10_02.vsh10_03"]] <- as.data.frame(cts[, c("Neuron_h10_02_bTx1", "Neuron_h10_03_bIGG")])

# counts_astro.neuron_mat <- lapply(counts_astro.neuron, as.matrix)

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(cts))
rownames(colData)<-NULL

colData <- colData %>%
  mutate(condition = ifelse(grepl('a5_02', samples), 'AD+Irisin', 'AD')) %>%
  column_to_rownames(var = 'samples')

# get more information from metadata

# sum.count <- counts_astro.neuron[["astro_h10_01.vsh10_02"]]


# perform DESeq2 --------
# Create DESeq2 object   
# YJ.Edit
# deseq.fn <- function(x){x <- DESeqDataSetFromMatrix(countData = x,
#                                                     colData = colData,
#                                                     design = ~ condition)}
# 
# as.n.deseq <- lapply(counts_astro.neuron, deseq.fn)

# F.count.data = sum.count[sum.count >5,]
# count.mat = F.count.data
# 
# conds = DataFrame(condition = factor(c("K","L")))
# dds = DESeqDataSetFromMatrix(F.count.data, conds, formula(~ condition))
# 
# dds = DESeq(dds)
# res = results(dds)

library(DESeq2)
library(tidyverse)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)



# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_AD.Irisin_vs_AD")
res

write.table(res, file='condition_AD.Irisin_vs_AD.tsv', quote=FALSE, sep='\t', col.names = TRUE)
# PseudoBulk end ----

seurat_Ast_A5.AB.ast <- subset(seurat_Ast_A5.AB, idents = c("Astrocyte"))
features <- c("RELA", "GFAP" , "BDNF", "MME", "NFKB1", "NFKB2")
VlnPlot(seurat_Ast_A5.AB, features = features)
help("VlnPlot")

# VlnPlot export ----
# tri.subset.CD8_WT <- subset(tri.subset.CD8, idents = "Immune_WT")
# write.table(tri.subset.CD8_WT@assays[["RNA"]]@data, file='tri.subset.CD8_WT_Gene_data_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)
# tri.subset.CD8_AD <- subset(tri.subset.CD8, idents = "Immune_AD")
# write.table(tri.subset.CD8_AD@assays[["RNA"]]@data, file='tri.subset.CD8_AD_Gene_data_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)


 


vln.seurat_Ast_A5.AB <- subset(seurat_integrated.3D.ast, idents = c( "a5_01_b",   "a5_02_b",   "a5_03_b"))
features <- c("RELA", "GFAP" , "BDNF", "MME", "NFKB1", "NFKB2")
VlnPlot(vln.seurat_Ast_A5.AB, features = features)

FeaturePlot(vln.seurat_Ast_A5.AB, features = features)

DotPlot(vln.seurat_Ast_A5.AB, features = features) + RotatedAxis()
# Single cell heatmap of feature expression
DoHeatmap(subset(vln.seurat_Ast_A5.AB, downsample = 100), features = features, size = 3)

# export as xlsx
library(xlsx)
for(i in 1:length(de.line.list_rna)){
  write.xlsx(de.line.list_rna[[i]], 
             file = "results/Ast_de_line_rna_results.2022.05.20.xlsx", 
             sheetName = names(de.line.list_rna[i]),
             append = T)}

# Neuron
seurat_integrated.3D.neu <- subset(seurat_integrated.selected.01.combined.name, idents = c("Neuron"))
# Neuron DE by line
Idents(seurat_integrated.3D.neu) <- "sample"

de.line.list_rna <- list()
de.line.list_rna[["a5_02_a vs. a5_01_a"]] <- de(seurat_integrated.3D.neu, 
                                                  id1 = "a5_02_a", 
                                                  id2 = "a5_01_a", 
                                                  assay = "RNA")
de.line.list_rna[["a5_03_a vs. a5_02_a"]] <- de(seurat_integrated.3D.neu, 
                                                  id1 = "a5_03_a", 
                                                  id2 = "a5_02_a", 
                                                  assay = "RNA")


# export as xlsx
library(xlsx)
for(i in 1:length(de.line.list_rna)){
  write.xlsx(de.line.list_rna[[i]], 
             file = "results/Neu_de_line_rna_results.2022.05.20.xlsx", 
             sheetName = names(de.line.list_rna[i]),
             append = T)}



Neu <- subset(seurat_integrated, idents = c("5", "7"))

de <- function(x, id1, id2, assay){
  x <- FindMarkers(x,
                   ident.1 = id1,
                   ident.2 = id2,
                   test.use = "wilcox",
                   assay = assay,
                   slot = "data",
                   min.pct = 0.05,
                   logfc.threshold = 0.05)
  return(x)
}

# Astrocyte DE by line
Idents(Ast) <- "sample"

de.line.list_rna <- list()
de.line.list_rna[["h10_02_b vs. h10_01_b"]] <- de(Ast, 
                                                  id1 = "h10_02_b", 
                                                  id2 = "h10_01_b", 
                                                  assay = "RNA")
de.line.list_rna[["h10_03_b vs. h10_02_b"]] <- de(Ast, 
                                                  id1 = "h10_03_b", 
                                                  id2 = "h10_02_b", 
                                                  assay = "RNA")


# export as xlsx
library(xlsx)
for(i in 1:length(de.line.list_rna)){
  write.xlsx(de.line.list_rna[[i]], 
             file = "results/Ast_de_line_rna_results.xlsx", 
             sheetName = names(de.line.list_rna[i]),
             append = T)}



DimPlot(seurat_integrated, pt.size = 1, label = TRUE, label.size = 7)
  
  
  
  
# 2022.05.20 ByJoseph end ----












































































































# comeback here ====
# Explore resolutions, there are different columns for each of the resolutions calculated.
seurat_integrated@meta.data %>% 
  View()


# _ Clusters by Resolution ----
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4.5) +
  ggtitle("Clustering at 0.2") + 
  theme(plot.title = element_text(hjust=0.5, face="bold"))
  
seurat_integrated$integrated_snn_res.0.2

# _ Clustering QC ----
Idents(seurat_integrated) <- "stim"
# Segregation of clusters by sample: cell distribution per cluster in each sample
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "integrated_snn_res.0.2")) %>%
  dplyr::count(ident, integrated_snn_res.0.2) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()


# Segregation of clusters by various sources of uninteresting variation
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
DefaultAssay(seurat_integrated) <- "RNA"
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',             # play around with these number (min.cutoff, max.cutoff)
            # max.cutoff = 'q99',
            label = TRUE,
            ncol = 3)


# . ----
# 10. Cell Type Identification ----
# 10-1. FeaturePlot----
DefaultAssay(seurat_integrated) <- "RNA"

fplot <- function(markers, ncol, min, max){FeaturePlot(seurat_integrated,
                                                       reduction = "umap",
                                                       features = markers,
                                                       label = TRUE,
                                                       col = c("grey", "blue"),
                                                       ncol = ncol,
                                                       min.cutoff = NA,
                                                       max.cutoff = NA)
}


# gene markers
m.astro <- c("AQP4", "UBE2C", "NUSAP1", "TOP2A", "PTPRZ1", "HOPX", "FAM107A")
# general astro: "AQP4"
# fetal astro: "UBE2C", "NUSAP1", "TOP2A"
# intermediate astro: "PTPRZ1", "HOPX", "FAM107A"
m.neurons <- c("ABAT", "GAD1", "KCNJ6", "TPH1", "DCX", "GABBR2", "SATB2", "FOXP2")
m.oligo <- c("FTH1", "PLP1", "FGF1", "MBP", "MOBP")
m.opcs <- c("PDGFRA", "VCAN", "CSPG4")
# m.mic <- c("HLA-DRB1")
markers.all <- c(m.astro, m.neurons, m.oligo, m.opcs)
select.markers <- c("AQP4", "AGT", "GABBR2", "KCNJ6")

fplot(m.astro, 3, "q10", "q90")
fplot(m.neurons, 3)
fplot(m.neurons, 3, min = "q25")
fplot(m.oligo, 3)
fplot(m.oligo, 5, min = "q25", "q90")
fplot(m.opcs, 2)



# 10-2. DotPlot ----
# color setting. 
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9, "YlOrRd")
my_color <- brewer.pal(9, "YlOrRd")
my_color <- c('darkgoldenrod1', 'orange', 'turquoise3', 'red')
my_color


View(seurat_integrated)
dot <- function(gene, title){
  DotPlot(seurat_integrated, features = gene, dot.scale = 6) +
    RotatedAxis() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5))  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}

dot(markers.all, "DotPlot_ALL_resolution 0.2")
dot(m.astro, "DotPlot_Astrocyte_resolution 0.2")
dot(m.neurons, "DotPlot_Neurons_resolution 0.2")
dot(m.oligo, "DotPlot_Oligodendrocyte_resolution 0.2")
dot(m.opcs, "DotPlot_OPCs_resolution 0.2")



# 10-3. Conserved Markers----
# To identify canonical cell type marker genes that are conserved across conditions
# because cluster 21 has only 1 cell in ctrl, so doing conserved markers only for cluster0 through 20. For cluster 21, findmarkers will be performed.
c.markers <- list()
for(i in c(0:20)){
  c.markers[[i+1]] <- FindConservedMarkers(seurat_integrated, ident.1 = i,
                                         grouping.var = "stim")
}

names(c.markers) <- paste0("c.markers_", 0:20)

# FindMarkers for cluster 21 (not grouped by stim)
cluster21.markers <- FindMarkers(seurat_integrated, ident.1 = 21)
head(cluster21.markers)
cluster21.markers %>% nrow()



# comeback here ====


# Extract top n convserved markers per cluster
library(tibble)
library(dplyr)
top_n <- function(x, n){x %>% 
    rownames_to_column(var = "gene") %>%
    mutate(avg_logFC = (AD_avg_logFC + CTRL_avg_logFC)/2) %>% 
    slice_max(avg_logFC, n = n) %>%
    select(gene, max_pval, minimump_p_val, avg_logFC)
}

c.markers_top30 <- lapply(c.markers, top_n, n=30)

# export
library(xlsx)
for (i in 1:length(c.markers_top30)){
  write.xlsx(c.markers_top30[[i]], 
             file="results/conserved.markers_1_20.xlsx", 
             sheetName=names(c.markers_top30[i]), 
             append=T,
             row.names = F)
}


# see what's top genes in cluster 21. play with the p-val_adj and how many genes to see
cluster21.markers %>% 
  filter(p_val_adj < 0.5) %>%
  slice_max(avg_logFC, n=nrow(cluster21.markers)) %>%
  View()



# 11. Re-clustering ----
seurat.select <- subset(seurat_integrated, idents = c("1", "3", "4", "21"))


# Re-do SCTransform before re-clustering (Normalize and CellCycleScoring is for the purpose of QC only)
# perform a log normalization
seurat.select <- NormalizeData(seurat.select, verbose = TRUE)
# cell cycle scoring for later (added to the metadata)
seurat.select <- CellCycleScoring(seurat.select, g2m.features=g2m_genes, s.features=s_genes)
# real data pre-processing done here. generating the SCT Assay
seurat.select <- SCTransform(seurat.select, vars.to.regress = c("percent.mt", "batch"))


# Reductions
# Run PCA
seurat.select <- RunPCA(object = seurat.select)


# Plot PCA
PCAPlot(seurat.select, split.by = "stim")
PCAPlot(seurat.select, split.by = "sample")
PCAPlot(seurat.select, split.by = "batch")


# Run UMAP
seurat.select <- RunUMAP(seurat.select, 
                         dims = 1:40,
                         reduction = "pca")


seurat.select <- RunTSNE(seurat.select, 
                         dims = 1:40,
                         reduction = "pca")


# Re-clustering
# Determine the K-nearest neighbor graph
seurat.select <- FindNeighbors(object = seurat.select, 
                               dims = 1:40)


# Determine the clusters for various resolutions                                
seurat.select <- FindClusters(object = seurat.select,
                              resolution = c(0.1, 0.2, 0.4))


# Explore resolutions, there are different columns for each of the resolutions calculated.
seurat.select@meta.data %>% 
  View()

View(seurat.select)
table(seurat.select$SCT_snn_res.0.1)
table(seurat.select$SCT_snn_res.0.2)
table(seurat.select$SCT_snn_res.0.4)

# _ Clusters by Resolution ----
Idents(seurat.select) <- "SCT_snn_res.0.1"


DimPlot(seurat.select,
        reduction = "umap",
        label = TRUE,
        label.size = 4.5)


#_ Clustering QC ----
meta.select <- seurat.select@meta.data

Idents(seurat.select) <- "stim"
Idents(seurat.select) <- "batch"
# Segregation of clusters by sample: cell distribution per cluster in each sample
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat.select, 
                     vars = c("ident", "SCT_snn_res.0.1")) %>%
  dplyr::count(ident, SCT_snn_res.0.1) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)






# UMAP of cells in each cluster by sample
Idents(seurat.select) <- "SCT_snn_res.0.1"
DimPlot(seurat.select, 
        label = TRUE, 
        split.by = "stim")

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat.select,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

DimPlot(seurat.select, 
        label = TRUE, 
        split.by = "batch")

# Segregation of clusters by various sources of uninteresting variation
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
DefaultAssay(seurat.select) <- "RNA"
fp.rna <- FeaturePlot(seurat.select, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',             # play around with these number (min.cutoff, max.cutoff)
            # max.cutoff = 'q99',
            label = TRUE)


DefaultAssay(seurat.select) <- "integrated"
fp.int <- FeaturePlot(seurat.select, 
                      reduction = "umap", 
                      features = metrics,
                      pt.size = 0.4, 
                      order = TRUE,
                      min.cutoff = 'q10',             # play around with these number (min.cutoff, max.cutoff)
                      # max.cutoff = 'q99',
                      label = TRUE)



fp.rna
fp.int

# .----
# 12. Cell type identification ----
# 12-1. FeaturePlot----
DefaultAssay(seurat.select) <- "RNA"

fplot.2nd <- function(markers, ncol, min, max){FeaturePlot(seurat.select,
                                                           reduction = "umap",
                                                           features = markers,
                                                           label = TRUE,
                                                           col = c("grey", "blue"),
                                                           ncol = ncol,
                                                           min.cutoff = NA,
                                                           max.cutoff = NA)
}


fplot.2nd(m.astro, 4)
fplot.2nd(m.neurons, 4)
fplot.2nd(select.markers, 4)


Idents(seurat.select) <- "SCT_snn_res.0.1"
Idents(seurat.select) <- "SCT_snn_res.0.2"
fplot.2nd(c(m.astro, m.neurons), 4)



# 12-2. DotPlot ----
Idents(seurat.select) <- "SCT_snn_res.0.1"
DotPlot(seurat.select, features = c(m.astro, m.neurons), dot.scale = 6) +
  RotatedAxis() +
  labs(title = "DotPlot_All_resolution 0.1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))


DotPlot(seurat.select, features = select.markers, dot.scale = 6) +
  RotatedAxis() +
  labs(title = "DotPlot_All_resolution 0.2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))


Idents(seurat.select) <- "SCT_snn_res.0.2"
DotPlot(seurat.select, features = c(m.astro, m.neurons), dot.scale = 6) +
  RotatedAxis() +
  labs(title = "DotPlot_All_resolution 0.2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))







# . ----
# . ----
# Analysis on h10 batch b data ----
seurat_b_h10 <- seurat.list[3:5]

# save this as a permanent file and load them
save(seurat_b_h10, file="14_01. B_H10/seurat objects/b_h10.RData")
# load("seurat objects/b_h10.RData")
# The analysis should be held in a separate project (14_01. B_H10)




# . ----
# . ----
# Monocle 3 ----
# 12/7, 12/29, & 2/8

# 0. Set up ----
# Installation
# Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

# Install a few Bioconductor dependencies that aren't automatically installed
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

# install monocle3 through the cole-trapnell-lab GitHub
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase', force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3', force = TRUE)

library(monocle3)


# 1. Creat CDS ----
# Using split_seurat (split by sample based on filtered Seurat objects)
# Not possible to convert Seurat objects to CDS because importCDS is not built in monocle3 yet. So need to build CDS using each element


# make a list of counts matrix
count.mat_list <- list()
names(count.mat_list) <- names(split_seurat)

get.matrix <- function(x){
  as.matrix(GetAssayData(x, slot = "counts"))
}

get.sparse <- function(x){
  as.sparse(GetAssayData(x, slot = "counts"))
}

count.mat_list <- lapply(split_seurat, get.sparse)


# Cell metadata
# extract the object by sample and take only metadata
cell.meta.list <- list()
for(i in 1:length(split_seurat)){
  cell.meta.list[[i]] <- split_seurat[[i]]@meta.data
}

names(cell.meta.list) <- paste0("cell.meta_", names(split_seurat))

unlist(lapply(cell.meta.list, nrow))


# Gene annotation
gene.anno_list <- list()
for(i in 1:length(split_seurat)){
  gene.anno_list[[i]] <- data.frame("gene_short_name" = rownames(split_seurat[[i]]))
  rownames(gene.anno_list[[i]]) <- gene.anno_list[[i]][, 1]
}


names(gene.anno_list) <- names(split_seurat)


# Assembly 
cds.list <- list()
for(i in 1:length(split_seurat)){
  cds.list[[i]] <- new_cell_data_set(count.mat_list[[i]],
                                     cell_metadata = cell.meta.list[[i]],
                                     gene_metadata = gene.anno_list[[i]])
}

names(cds.list) <- names(split_seurat)


# Combining CDS objects
combined_cds <- combine_cds(cds.list)

# Combining H10_b CDS objects
combined_cds_b <- combine_cds(cds.list[3:5])

# normalization/pre-process + batch effect removal
pre_batch_align <- function(x){
  x <- preprocess_cds(x, num_dim = 100)
  x <- align_cds(x, alignment_group = "batch")
}



# normalization/pre-process + align with batch effect removal & mt content
pre_align <- function(x){
  x <- preprocess_cds(x, num_dim = 100)
  x <- align_cds(x, alignment_group = "batch", residual_model_formula_str = "~percent.mt")
}

combined_cds._h10b_final <- pre_align(combined_cds_b)


# rename the object for convenience
cds_h10b <- combined_cds._h10b_final

# Check explained variance
plot_pc_variance_explained(cds_h10b)

# Reduce dimensionality
cds_h10b <- reduce_dimension(cds_h10b)
plot_cells(cds_h10b)
plot_cells(cds_h10b, genes=c("AQP4", "GFAP"))
plot_cells(cds_h10b, genes = "STAT1")
plot_cells(cds_h10b, genes = "AQP4")
plot_cells(cds_h10b, color_cells_by = "batch")

# Group cells into clusters
cds_h10b <- cluster_cells(cds_h10b, resolution = 1e-2)
plot_cells(cds_h10b)

cds_h10b <- cluster_cells(cds_h10b, resolution = 1e-3)
plot_cells(cds_h10b)
plot_cells(cds_h10b, genes = "AQP4")

plot_cells(cds_h10b, color_cells_by="partition", group_cells_by="partition")

# Chosse subcluster
cds_h10b_Ast <- choose_cells(cds_h10b)

# save permenant
saveRDS(cds_h10b_Ast, file = "cds objects/cds_h10b_Ast.rds")

cds_h10b_Ast_pr_graph_test_res <- graph_test(cds_h10b_Ast, neighbor_graph="knn", cores=8)
cds_h10b_Ast_pr_deg_ids <- row.names(subset(cds_h10b_Ast_pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

cds_h10b_Ast_gene_module_df <- find_gene_modules(cds_h10b_Ast[cds_h10b_Ast_pr_deg_ids,], resolution=1e-3)

plot_cells(cds_h10b_Ast, genes=cds_h10b_Ast_gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

cds_h10b_Ast = cluster_cells(cds_h10b_Ast, resolution=1e-3)
plot_cells(cds_h10b_Ast, color_cells_by="cluster")

#Learn the trajectory graph ====
cds_h10b_Ast <- learn_graph(cds_h10b_Ast)
plot_cells(cds_h10b_Ast,
           color_cells_by = "drug",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds_h10b_Ast <- order_cells(cds_h10b_Ast)
plot_cells(cds_h10b_Ast,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


ciliated_genes <- c( "RELA"
                    )


cds_h10b_Ast_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

View(cds_h10b_Ast_subset)
table(cds_h10b_Ast_subset@colData$drug)
cds_h10b_Ast_subset@colData$drug <- factor(cds_h10b_Ast_subset@colData$drug, levels = c("Ctrl", "Tx1", "IGG"))

# save permanent
saveRDS(cds_h10b_Ast_subset, file = "cds objects/cds_h10b_Ast_subset.rds")
# cds_h10b_Ast_subset <- readRDS(file = "cds objects/cds_h10b_Ast_subset.rds")

gene_fits <- fit_models(cds_h10b_Ast_subset, model_formula_str = "~drug")
fit_coefs <- coefficient_table(gene_fits)

library(dplyr)
emb_time_terms <- fit_coefs %>% filter(term == "drug")
emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

# plot_genes_violin(cds_h10b_Ast_subset, group_cells_by="drug", ncol=2)
plot_genes_violin(cds_h10b_Ast_subset, group_cells_by="drug", ncol=3) +
  theme(axis.text.x=element_text(angle=45, hjust=1))



# ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
# pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))


# **************** DE cds_h10b_Ast**************** ----
# * DE, Regression model ----
# Ctrl, _02: Tx1, 03_IGG
# Re-leveling of metadata drug
levels(factor(cds_h10b_Ast@colData$drug))
cds_h10b_Ast@colData$drug <- factor(cds_h10b_Ast@colData$drug, levels = c("Ctrl", "Tx1", "IGG"))

imp_genes <- c("GFAP", "S100B", 
               "C3", "C3AR", 
               "STAT3", "NFKB1", "IL6",
               "MME", "MAPK3",
               "ITGAV", "ITGB5", "RELA",
               "MT-CO1", "PTGS2")


imp_genes <- c("GFAP", "S100B")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=4)

imp_genes1 <- c("C3", "C3AR")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)

imp_genes <- c("STAT3", "NFKB1", "IL6")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)

imp_genes <- c("MME", "MAPK3")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)

imp_genes <- c("ITGAV", "ITGB5", "RELA")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)

imp_genes <- c("MT-CO1", "PTGS2")
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)
               
cds_h10b_Ast <- cds[rowData(cds)$gene_short_name %in% imp_genes,]

# or you can pull the rows by using rownames if you don't have gene_short_name in the following way
# cds_subset <- cds.noNaive[rownames(cds.noNaive) %in% imp_genes, ]

gene_fits <- fit_models(cds_h10b_Ast, model_formula_str = "~drug")

# View(gene_fits)

# see which of the genes have line-dependent expression
fit_coefs <- coefficient_table(gene_fits)
# View(fit_coefs)

# filtering out intercept terms (bc generally don't care about the intercept term, beta0)
line_terms <- fit_coefs %>% filter(term != "(Intercept)")
# View(line_terms)

# filtering based on q-val (you can edit: q_val <=0.05) 
line_terms_filtered <- line_terms %>% 
  filter (q_value < 0.05) %>% 
  select(gene_short_name, term, q_value, estimate)
# View(line_terms_filtered)

# Violin ====
plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2)

plot_genes_violin(cds_h10b_Ast, group_cells_by="drug", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# Pseudoplots - Finding genes that change as a function of pseudotime ====
plot_cells(cds_h10b_Ast,
           color_cells_by = "drug",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ciliated_cds_pr_test_res <- graph_test(cds_h10b_Ast, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

#Pseudotime
cds_h10b_Ast <- order_cells(cds_h10b_Ast)
AFD_genes <- c("GFAP", "S100B", 
               "C3", "C3AR", 
               "STAT3", "NFKB1", "IL6",
               "MME", "MAPK3",
               "ITGAV", "ITGB5", "RELA",
               "MT-CO1", "PTGS2")
AFD_lineage_cds <- cds_h10b_Ast[rowData(cds_h10b_Ast)$gene_short_name %in% AFD_genes, ]

AFD_lineage_cds <- cds_h10b_Ast[rowData(cds_h10b_Ast)$gene_short_name %in% AFD_genes,
                            colData(cds_h10b_Ast)$drug %in% c("Ctrl", "Tx1", "IGG")]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="drug",
                         min_expr=0.1)

plot_genes_in_pseudotime(AFD_lineage_cds,
                         
                         min_expr=0.001)






# **************** DE cds_h10b**************** ----
# * DE, Regression model ----
# Ctrl, _02: Tx1, 03_IGG

imp_genes <- c("GFAP", "S100B", 
               "C3", "C3AR", 
               "STAT3", "NFKB1", "IL6",
               "MME", "MAPK3",
               "ITGAV", "ITGB5", "RELA",
               "MT-CO1", "PTGS2")

cds_h10b <- cds[rowData(cds)$gene_short_name %in% imp_genes,]

# or you can pull the rows by using rownames if you don't have gene_short_name in the following way
# cds_subset <- cds.noNaive[rownames(cds.noNaive) %in% imp_genes, ]

gene_fits <- fit_models(cds_h10b, model_formula_str = "~drug")
gene_fits <- fit_models(cds_subset, model_formula_str = "~drug + line")
View(gene_fits)

# see which of the genes have line-dependent expression
fit_coefs <- coefficient_table(gene_fits)
View(fit_coefs)

# filtering out intercept terms (bc generally don't care about the intercept term, beta0)
line_terms <- fit_coefs %>% filter(term != "(Intercept)")
View(line_terms)

# filtering based on q-val (you can edit: q_val <=0.05) 
line_terms_filtered <- line_terms %>% 
  filter (q_value < 0.05) %>% 
  select(gene_short_name, term, q_value, estimate)

View(line_terms_filtered)

# Violin ====
library(ggplot2)
plot_genes_violin(cds_h10b, group_cells_by="drug", ncol=2)
plot_genes_violin(cds_h10b, group_cells_by="drug", ncol=4) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


# Pseudoplots - Finding genes that change as a function of pseudotime ====
plot_cells(cds_h10b,
           color_cells_by = "drug",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ciliated_cds_pr_test_res <- graph_test(cds_h10b, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

#Pseudotime
cds_h10b <- order_cells(cds_h10b)
AFD_genes <- c("GFAP",
               "STAT3",
               "NFKB1",
               "ITGAV",
               "MT-CO1",
               "ITGB5")
AFD_lineage_cds <- cds_h10b[rowData(cds_h10b)$gene_short_name %in% AFD_genes, ]

AFD_lineage_cds <- cds_h10b[rowData(cds_h10b)$gene_short_name %in% AFD_genes,
                               colData(cds_h10b)$drug %in% c("Ctrl", "Tx1", "IGG")]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="drug",
                         min_expr=0.5)

