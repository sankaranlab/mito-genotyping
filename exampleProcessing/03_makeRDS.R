library(data.table)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tools)
library(Matrix)
library(GenomicRanges)
library(seqinr)

#' Import data from plain text files explicitly point
#' to the requisite input files and render a RangedSummarized
#' Experiment.
#'
#' \code{importMito.explicit} takes a sparse matrix file of position,
#' sample, and count (non-zero) in some order as variant calls
#' (one file per allele) as well as
#' a separate file of the position and sample coverage and
#' produces a MultiAssayExperiment object that serves
#' as the backbone of the R interface 
#'
#' The filepath of a plain text or gzipped file
#' that contains data in a long matrix format of the position,
#' allele, sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' analyses
#'
#' @param Afile Contains the position, sample, count of A alleles
#' @param Cfile Contains the position, sample, count of C alleles
#' @param Gfile Contains the position, sample, count of G alleles
#' @param Tfile Contains the position, sample, count of T alleles
#'
#' @param coverageFile The filepath of a plain text or gzipped file
#' that contains data in a sparse matrix format of the position,
#' sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' analyses.
#'
#'
#' @param referenceAlleleFile The filepath of a plain text fasta
#' file of the mtDNA
#'
#'
#' @return An initialized object that is an S4 class
#' MultiAssayExperiment that includes variant
#' count and coverage as assay slots.
#'
#' @import Matrix
#' @importFrom data.table fread dcast.data.table
#' @importFrom tools file_ext
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import MultiAssayExperiment
importMito.explicit <- function(Afile, Cfile, Gfile, Tfile,
                                coverageFile, referenceAlleleFile){
  
  variantFiles <- list(Afile, Cfile, Gfile, Tfile)
  
  # Set up downstream processing including robust ordering
  # The coverage file could have slightly more variants /
  # individual samples depending on the calls, so base it
  # of of them
  importDT <- function(file){
    if(tools::file_ext(file) == "gz"){
      cov <- data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      cov <- data.table::fread(paste0(file), stringsAsFactors = TRUE)
    } else{
      stop("Provide a valid file format for the  file (.gz, .txt, .csv, or .tsv)")
    }
  }
  
  cov <- importDT(coverageFile)
  
  samplesOrder <- levels(cov[[2]])
  maxpos <- max(cov[[1]])
  maxsamples <- length(samplesOrder)
  
  # make coverage a sparse matrix
  covmat <- Matrix::sparseMatrix(
    i = cov[[1]],
    j = as.numeric(cov[[2]]),
    x = cov[[3]]
  )
  remove(cov)
  
  # Import Counts and BAQ
  importSMs <- function(file){
    # fread the individual variant calls in
    if(tools::file_ext(file) == "gz"){
      dt <-  data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      dt <-  data.table::fread(paste0(file), stringsAsFactors = TRUE)
    } else{
      stop("Provide a valid file format for the variant call file (.gz, .txt, .csv, or .tsv)")
    }
    
    dt$sample <- factor(dt$sample, levels = samplesOrder)
    
    counts <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[3]],0)
    )
    
    BAQ <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[4]],0)
    )
    remove(dt)
    return(list("counts" = counts, "BAQ" = BAQ))
  }
  
  ACGT <- lapply(variantFiles, importSMs)
  names(ACGT) <- c("A", "C", "G", "T")
  
  # Import fasta
  reffasta <- seqinr::read.fasta(fasta)
  mitoChr = attr(reffasta, "name")
  ref <- data.frame(
    V1 = 1:maxpos,
    V2 = toupper(unname(reffasta)[[1]])[1:maxpos]
  )
  
  
  # Make a long matrix of BAQ and Counts for non-reference alleles
  whichA <- which(ref[["V2"]][1:maxpos] != "A")
  whichC <- which(ref[["V2"]][1:maxpos] != "C")
  whichG <- which(ref[["V2"]][1:maxpos] != "G")
  whichT <- which(ref[["V2"]][1:maxpos] != "T")
  
  longBAQ <- rbind(
    ACGT[["A"]][["BAQ"]][whichA,],
    ACGT[["C"]][["BAQ"]][whichC,],
    ACGT[["G"]][["BAQ"]][whichG,],
    ACGT[["T"]][["BAQ"]][whichT,]
  )
  
  longCounts <- rbind(
    ACGT[["A"]][["counts"]][whichA,],
    ACGT[["C"]][["counts"]][whichC,],
    ACGT[["G"]][["counts"]][whichG,],
    ACGT[["T"]][["counts"]][whichT,]
  )
  letterz <- c(rep("A", length(whichA)), rep("C", length(whichC)), rep("G", length(whichG)), rep("T", length(whichT)))
  remove(ACGT)
  
  # Create colData
  depth <- Matrix::colMeans(covmat)
  sdf <- data.frame(sample = samplesOrder, depth)
  rownames(sdf) <- samplesOrder
  colnames(sdf) <- c("sample", "depth")
  
  # Make row Ranges for each object
  row_g_cov <- GenomicRanges::GRanges(seqnames = mitoChr,
                                      IRanges::IRanges(1:maxpos, width = 1))
  GenomicRanges::mcols(row_g_cov) <- data.frame(refAllele = ref[[2]][1:maxpos])
  
  row_g_allele <- GenomicRanges::GRanges(seqnames = mitoChr,
                                         IRanges::IRanges(1:maxpos, width = 1))[c(whichA, whichC, whichG, whichT)]
  
  GenomicRanges::mcols(row_g_allele) <- data.frame(refAllele = (ref[[2]][1:maxpos])[c(whichA, whichC, whichG, whichT)],
                                                   altAllele = letterz)
  
  # Make summarized experiments and
  coverage <- SummarizedExperiment::SummarizedExperiment(
    assays = list("coverage" = covmat),
    colData = S4Vectors::DataFrame(sdf),
    rowData = row_g_cov
  )
  
  alleles <- SummarizedExperiment::SummarizedExperiment(
    assays = list("BAQ" = longBAQ, "counts" = longCounts),
    colData = S4Vectors::DataFrame(sdf),
    rowData = row_g_allele
  )
  
  remove(longBAQ)
  remove(longCounts)
  
  MAE <- MultiAssayExperiment::MultiAssayExperiment(
    list("alleles" = alleles, "coverage" = coverage),
    colData = S4Vectors::DataFrame(sdf)
  )
  return(MAE)
}


#' Import data from folder output of preprocessing.
#'
#' \code{importMito} takes a path to a folder that contains
#' all the requisitive input files for setting up an analysis.
#'
#' @param folder Filepath to folder 
#'
#' @param ... Additional parameters to pass to the
#' importMito.explict function
#'
#' @return An initialized  object that is an S4 class
#' MultiAssayExperiment that includes variant
#' frequency and coverage as assay slots.
#'
#' @seealso importMito.explicit
importMito.folder <- function(folder, fasta){
  
  files <- list.files(folder, full.names = TRUE)
  
  checkGrep <- function(hit){
    if(length(hit) != 1){
      stop("Improper folder specification; file missing / extra file present. See documentation")
    } else {
      return(hit)
    }
  }
  
  # Set up file paths
  Afile <- files[checkGrep(grep(".A.txt.gz", files))]
  Cfile <- files[checkGrep(grep(".C.txt.gz", files))]
  Gfile <- files[checkGrep(grep(".G.txt.gz", files))]
  Tfile <- files[checkGrep(grep(".T.txt.gz", files))]
  coverageFile <- files[checkGrep(grep(".coverage.txt.gz", files))]
  
  referenceAlleleFile <- fasta
  
  
  MAE <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
                             coverageFile, referenceAlleleFile)
  return(MAE)
}

#folder = "processed_data"
#fasta = "mito_fastas/hg19.fasta"

args <- commandArgs(trailingOnly = FALSE)

# Parse parameters working backwards since I am keeping all of them
# to be able to source the knee call Rscript
nn <- length(args)
folder <- args[nn-1]
fasta <- args[nn]

MAE <- importMito.folder(folder, fasta)
saveRDS(MAE, paste0(folder, "/processed.MAE_mito.rds"))


