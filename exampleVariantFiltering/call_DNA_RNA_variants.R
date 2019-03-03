library(SummarizedExperiment)
library(ggplot2)
library(BuenColors)
library(MultiAssayExperiment)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(plotly)

"%ni%" <- Negate("%in%")

raw <- readRDS("../mito_data/atac_rna/pop3_ATAC-RNA.rds")

covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
baq <- assays(allSE)[["BAQ"]]
rm(raw)
covCell <-  Matrix::colMeans(assays(covSE)[["coverage"]])

# Allele Frequency
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]

master_df <- data.frame(cbind(data.matrix(af), data.matrix(cov), data.matrix(baq)))
colnames(master_df) <- c("AF_ATAC", "AF_RNA", "Cov_ATAC", "Cov_RNA", "ATAC_BQ", "RNA_BQ")
master_df$position <- data.frame(rowRanges(allSE))[,c(2)]
master_df$altAllele <- data.frame(rowRanges(allSE))[,c(7)]
master_df$refAllele <- data.frame(rowRanges(allSE))[,c(6)]

if(FALSE){
  
  h1 <- hist(master_df$ATAC_BQ, breaks=150, plot=FALSE)
  h2 <- hist(master_df$RNA_BQ, breaks=150, plot=FALSE)
  top <- max(h1$counts, h2$counts)
  k <- kde2d(master_df$ATAC_BQ, master_df$RNA_BQ, n=150)
  
  rf <- colorRampPalette(jdb_palette("brewer_spectra"))
  r <- rf(32)
  
  # margins
  pdf("../output/bivariate_BQ.pdf", height = 2.5, width = 2.5)
  oldpar <- par()
  par(mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
  image(k, col=r) #plot the image
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='black')
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='black', horiz=TRUE)
  dev.off()
  
}

boo <- (master_df$AF_ATAC > 0.001 | master_df$AF_RNA > 0.001) &
  (master_df$AF_ATAC < 0.99 | master_df$AF_RNA < 0.99) &
  (master_df$ATAC_BQ > 25) & (master_df$RNA_BQ > 25) &
  master_df$Cov_RNA > 10000 & master_df$Cov_ATAC > 1000

master_df$keepInAnalysis <- ifelse((sqrt(master_df$AF_ATAC) > sqrt(0.003)) &
                                     (master_df$AF_ATAC < 0.5 & master_df$AF_RNA < 0.5) &
                                     (master_df$ATAC_BQ > 20) & (master_df$RNA_BQ > 20) &
                                     master_df$Cov_RNA > 10000 & master_df$Cov_ATAC > 1000 , "Kept", "Filtered")

mmdf <- master_df[master_df$keepInAnalysis == "Kept",]
#mmdf$percent_difference <- (mmdf$AF_RNA - mmdf$AF_ATAC)/((mmdf$AF_RNA + mmdf$AF_ATAC)/2)*100
write.table(mmdf[order(mmdf$AF_ATAC, decreasing = FALSE),], file = "../output/allele_frequences_extended.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

if(FALSE){
  write.table(data.frame(cells = rownames(mmdf)), file = '../output/calledVariants.tsv', sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
}

p1 <- ggplot(arrange(master_df[boo, ], AF_ATAC), aes(x = sqrt(AF_ATAC), y = sqrt(AF_RNA), color = as.factor(keepInAnalysis), group = position))+
  geom_point(size = 0.5) + labs(x = "sqrt Allele Frequency - ATAC", y = "sqrt Allele Frequency - RNA", color = "Variant" ) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + pretty_plot(fontsize = 6) +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(legend.position = "none") + L_border()

ggplotly(p1)
cowplot::ggsave(p1, file = "../output/ATAC-RNA_variantCall.pdf", width = 1.7, height = 1.7)

pNoSR <- ggplot(arrange(master_df[boo, ], AF_ATAC), aes(x = (AF_ATAC), y = (AF_RNA), color = as.factor(keepInAnalysis), group = position))+
  geom_point(size = 0.5) + labs(x = "Allele Frequency - ATAC", y = "Allele Frequency - RNA", color = "Variant" ) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + pretty_plot(fontsize = 6) +
  scale_color_manual(values = c("black", "firebrick")) +
  theme(legend.position = "none") + L_border()

cowplot::ggsave(cowplot::plot_grid(p1, pNoSR, nrow = 1), file = "../output/ATAC-RNA_variantCall-R2R.pdf", width = 4, height = 2)
