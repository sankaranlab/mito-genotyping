library(SummarizedExperiment)
library(BuenColors)
library(MultiAssayExperiment)
library(gplots)
library(BioQC)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(plotly)
library(RColorBrewer)


"%ni%" <- Negate("%in%")

# Import Caleb's genotype calls
raw <- readRDS("../data/mgatk-data/P0322-mgatk/P0322-mgatk.rds")

# Import Lauren's clonal definitions; make a vector for quick lookup
meta_df <- data.frame(fread("../output/clones/P0322.txt"))
clonevec <- as.character(meta_df$Clonotype.ID); names(clonevec) <- as.character(meta_df$Cell.Name)

# Collapse multi assay experiments into SummarizedExperiments
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
#baq <- assays(allSE)[["BAQ"]] # don't need this for now.

# Pull out the allele frequencies and the coverages per position / allele
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], data.frame(rowRanges(allSE))[,c(6)], ">", data.frame(rowRanges(allSE))[,c(7)])

# Compute per-cell and per-variant properties
covVARmean <- Matrix::rowMeans(assays(covSE)[["coverage"]][start(rowRanges(allSE)),])
covVARmedian <- matrixStats::rowMedians(data.matrix(assays(covSE)[["coverage"]][start(rowRanges(allSE)),]))
covCell <-  Matrix::colMeans(assays(covSE)[["coverage"]])

# Filter variants and cells that don't have great coverages
minCellCoverage <- 50
af2 <- af[covVARmean> 20  & covVARmedian > 20 , covCell > minCellCoverage]
cov2 <- cov[covVARmean> 20 & covVARmedian > 20, covCell > minCellCoverage]
dim(af2)

# Filter for the minor allele
af3 <- af2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]
cov3 <- cov2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]
dim(af3)

# Adjust column names
colnames(af3) <- gsub(".mito", "", colnames(af3))
colnames(cov3) <- gsub(".mito", "", colnames(cov3))

# Remove samples where we don't have a clone ID... shouldn't be any
af3 <- af3[,colnames(af3) %in% names(clonevec)]
cov3 <- cov3[,colnames(cov3) %in% names(clonevec)]
dim(af3)

# Create a vector per cell with clone ID
clones <- clonevec[colnames(af3)]

# Identify all clones with at least 5 cells
data.frame(clones) %>% group_by(clones) %>%
  summarise(count = n()) %>% filter(count >=5) %>%
  pull(clones) %>% as.character() -> search_clones

# Loop over these clones to identify variants
lapply(search_clones, function(clone){
  
  # Do a Mann-Whitney test
  out <- sort(wmwTest(data.matrix(t(af3)), which(clones == clone),
                      valType = "abs.log10p.two.sided"), decreasing = TRUE)
  
  # Retain if reasonable p-value-- will threshold harder, later
  keep <- out[out> 4]
  df <- data.frame(Variant = names(keep), log10p = unname(keep), colony = clone)
  df
})  %>% rbindlist() %>% as.data.frame() -> variantScoresClones

# Verify that true positive variants are there
variantScoresClones %>% filter(Variant %in% c("3515T>C","13412T>C","12008G>A"))

# Try different threshold to see what is optimal
variantScoresClones %>% filter(log10p > 10) %>% dim()
keepVariants <- unique(variantScoresClones[variantScoresClones$log10p > 10, "Variant"])
length(keepVariants)

# Only keep variants that match this
boo <- rownames(af3) %in% keepVariants
af3 <- data.matrix(af3[boo,])
dim(af3)

# Threshold allele frequencies for a plot -- higher heteroplasmy values are definitely interesting
# but we will look at those case-by-case in box plots like in the existing manuscript
af4 <- af3
af4[af4 > .2] <- .2

# Make an ordered meta-data file
meta_df %>% filter(Cell.Name %in% colnames(af3)) %>% 
  mutate(clone_color = ifelse(Component.Size >= 5, Clonotype.ID, "small_clone")) %>%
  arrange(clone_color) -> annodf2
unique_clone_colors <- unique(as.character(annodf2$clone_color))

# Make a random vector of colors for annotation
colorsUniverse <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- sample(colorsUniverse, length(unique_clone_colors))
names(cols) <- unique_clone_colors

# That's it! Now make a heatmap
ha1 <- HeatmapAnnotation(df = annodf2[,c("clone_color"), drop = FALSE],
                         col = list(clone_color = cols)
)

pdf(file="../output/plots/heatmaps/P0322-heatmaps/P0322.heatmap.pdf", height = 9, width = 15)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(sqrt(af4[,as.character(annodf2$Cell.Name)]),
        col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha1,
        name = "AF")
dev.off()

# New heatmap but with actual allele frequencies
pdf(file="../output/plots/heatmaps/P0322.small.heatmap.pdf", height = 5, width = 15)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap((af3[c("3515T>C","13412T>C","12008G>A"),as.character(annodf2$Cell.Name)]),
        col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha1,
        name = "AF")
dev.off()
