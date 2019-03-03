library(SummarizedExperiment)
library(irlba)
library(BuenColors)
library(MultiAssayExperiment)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# This script reproduces Figure 3D in Ludwig et al

"%ni%" <- Negate("%in%")

# Import data
raw <- readRDS("tf1_barcode.rds")
meta <- read.table("full_sumstats_TF1barcode.tsv", header = TRUE)
meta <- meta[meta$group %ni% c("unassigned") & !is.na(meta$group),]
variants <- read.table("filtered_variants.tsv", header = TRUE, stringsAsFactors = FALSE)[,1]

# Unpack MAE and visualize allele frequency
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])
colnames(af) <- gsub(".mt", "", colnames(af))
af3 <- data.matrix(af[variants,as.character(meta$Cell)])

# Annodf 
annodf <- annodf2 <- meta

# Supervised detection of variants
groupsPossible <- unique(names(table(meta$group))[table(meta$group) > 3])
groups <- as.character(meta$group[meta$mtCoverage >100])
groups <- ifelse(is.na(groups), "noBarcode", groups)

# Make annotation
cols <- jdb_palette("lawhoops")[c(1, 3, 5, 7, 9, 11, 2, 13, 15, 17, 19)]; names(cols) <- unique(groups)

ha1 <- HeatmapAnnotation(df = annodf2[,c("group"), drop = FALSE],
                         col = list(group = cols)
)
af3[af3 > .2] <- .2
af3[af3 < 0.01] <- 0

pdf(file="calledHeatmap_cluster.pdf", width = 7, height = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(sqrt(af3), col=as.character(c("white", jdb_palette("brewer_red",type="continuous"))),
        cluster_rows = FALSE, cluster_columns = TRUE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize =4),
        top_annotation = ha1,
        name = "Allele\nFrequency")
dev.off()
