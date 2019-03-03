library(SummarizedExperiment)
library(irlba)
library(Rtsne)
library(ggplot2)
library(BuenColors)
library(MultiAssayExperiment)
library(ape)
library(gplots)
library(ggtree)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
source("00_distance.R")

"%ni%" <- Negate("%in%")

raw <- readRDS("../input/mito/atac_TF1.rds")
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]

baq <- assays(allSE)[["BAQ"]]
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)

if(FALSE){
  pp <- qplot(1:69, sort(colMeans(assays(covSE)[["coverage"]][start(rowRanges(allSE)),])), size = I(0.2)) +
    labs(x = "", y = "") + pretty_plot(fontsize = 4) + L_border + theme(legend.position = "none")
  cowplot::ggsave(pp, filename = "../output/final_plots/sampleCoverage.pdf", height = 1, width = 1)
}

baqdf <- data.frame(as.data.frame(rowRanges(allSE)), data.matrix(baq))

# Get non-zero row means
baqNonZero <- Matrix::rowSums(baq)/(Matrix::rowSums(baq > 0) + 0.001)
library(mixtools)
mixmdl <- normalmixEM(sort(baqNonZero), k = 3)

BAQcutoff <- min(mixmdl$x[mixmdl$posterior[,3] > 0.99])

pdf("../output/final_plots/BAQ_mixtureModel.pdf", height = 4, width = 4)
plot(mixmdl,which=2, breaks = 100)
lines(density(baqNonZero), lty=2, lwd=1)
abline(v = BAQcutoff, col="blue", lwd=1, lty=2)
mixmdl$mu
mixmdl$sigma
dev.off()

rm(baq)
rm(raw)
odf <- data.frame(rowRanges(allSE))
odf$rowMaxAF <- rowMax(data.matrix(af))
odf$BAQnonZero <- baqNonZero
odf$AF <- rowSums(assays(allSE)[["counts"]])/rowSums((assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001))
odf$coverage <- rowSums((assays(covSE)[["coverage"]][start(rowRanges(allSE)),]))

fiftydf <- odf %>% group_by(refAllele, altAllele) %>%
  summarise(fifty=quantile(BAQnonZero, probs=0.50), ten = quantile(BAQnonZero, probs=0.10)) 

odf_keep <- odf %>% group_by(refAllele, altAllele) %>%
  mutate(fifty=quantile(BAQnonZero, probs=0.50)) %>% ungroup() %>% arrange(start)
odf_keep$keep <- (odf_keep$fifty < odf_keep$BAQnonZero) & (BAQcutoff < odf_keep$BAQnonZero)
maxPos <- max(odf_keep$start)
odf_keep$keepAF <- (odf_keep$keep | (!odf_keep$keep & odf_keep$rowMaxAF < 0.002)) & (maxPos - 10 > odf_keep$start & odf_keep$start > 10)

valsVec <- tapply(odf_keep$keepAF, list(odf_keep$start), sum)

odf_keep[odf_keep$start == 2956,]
odf_keep[odf_keep$start == 1316,]
odf_keep[odf_keep$start == 1672,]
odf_keep[odf_keep$start == 2619,]

if(FALSE){
  library(ggbeeswarm)
  odf_keep$sqrtAF <- sqrt(odf_keep$AF)
  df <- odf_keep[odf_keep$BAQnonZero > 0,]
  df$alpha <- ifelse(df$keepAF, 1, 0.3)
  p1 <- ggplot(df[order(df$sqrtAF),], aes(x = altAllele, y = BAQnonZero, color = sqrtAF)) +
    geom_quasirandom(dodge.width=0.8, size = 2) +
    scale_color_gradientn(colors = jdb_palette("brewer_spectra")) +
    facet_wrap(~refAllele, strip.position = "bottom", nrow = 1) +theme_bw() +
    geom_point(data=fiftydf, aes(x = altAllele, y=fifty), color="black", shape = 95, size = 8, inherit.aes = FALSE) +
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(),
          strip.placement = "outside") + geom_hline(yintercept = BAQcutoff, color = "blue", linetype=2) +
    labs(x = "Reference Allele", y = "Mean BQ",  color = "Square Root\nAllele Frequency") +
    theme(legend.position = "bottom")
  
  cowplot::ggsave(p1, file = "../output/final_plots/bqAllele_plots.pdf", width = 7, height = 3)
}


# Allele Frequency
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]

# Filter using the Gini Index for generation 1
#sdf <- data.frame(stringr::str_split_fixed(colnames(af), pattern = "_", n = 5))
#af_g1 <- data.matrix(af[,sdf$X4 == "G001"])
#ginivals <- giniRows(af_g1)
#qplot(ginivals) + pretty_plot() + labs(x = "Gini Index", y = "count")

# Filter out variants that aren't get
# Don't go lower than 0.7 based on spot checking
minAF <- 0.025
maxAF <- 0.3
allele_freq <- rowMeans(data.matrix(af))
af3 <- af[as.character(start(rowRanges(allSE))) %in% names(valsVec[valsVec == 3]) & rowMaxs(data.matrix(af)) > minAF & allele_freq < maxAF, ]
cov3 <- cov[as.character(start(rowRanges(allSE))) %in% names(valsVec[valsVec == 3]) & rowMaxs(data.matrix(af)) > minAF & allele_freq < maxAF,]
dim(af3)

whichDF <- data.frame(pos = sapply(strsplit(rownames(af3), "_"), function(x) x[1]))
write.table(whichDF, file = "../output/whitelist_CL.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(data.frame(x = rownames(af3)), file = "../output/whitelist_CL_wLetter.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


df <- data.frame(stringr::str_split_fixed(colnames(af3), pattern = "_", n = 5))

tf1_color_maps2 <- tf1_color_maps[c(1,2,3, 13:16, 17, 10, 8, 7, 9, 6, 5 ,4, 11, 12)]

annodf <- df[,c("X3"), drop = FALSE]
colnames(annodf) <- "clone"
colnames(af3) <- df$X1
ha1 <- HeatmapAnnotation(df = annodf,
                         col = list(clone = tf1_color_maps2)
)

# Makea plot allele frequency graph
AT <- sapply(strsplit(rownames(af3), "_"), function(x) x[2]) %in% c("A","T")
paf <- af3
write.table(data.frame(rownames(paf)), file = "../output/whiteList_alleles.txt", sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

paf[paf > 0.2] <- 0.2
paf[paf < 0.0025] <- 0
pdf(file="../output/final_plots/tf1_heatmap_BQ.pdf", width = 5.5, height = 4)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(sqrt(data.matrix(paf)), col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 4.5),
        column_names_gp = gpar(fontsize = 4),
        top_annotation = ha1,
        name = "sqrt\nAllele\nFrequency")
dev.off()

#dim(data.matrix(af3))

if(FALSE){
  saveRDS(data.matrix(af3), file = "../output/alleleFreqs_bulkATAC.rds")
  saveRDS(data.matrix(cov3), file = "../output/alleleFreqs_coverages_bulkATAC.rds")
  twoMat <- distMito(af3, cov3)
  dist <- twoMat[[1]]
  saveRDS(dist, file = "../output/sample_distance.rds")
}
