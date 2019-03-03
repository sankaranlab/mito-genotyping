# Mitochondria genotyping
Repository for scripts for Ludwig et al. mtDNA genotyping as performed [in this manuscript](https://www.cell.com/cell/fulltext/S0092-8674(19)30055-8).

### Overview
Rather than performing variant calling using an established tool like GATK, 
this analytical framework instead 1) calculates all possible mutations' heteroplasmy
and 2) uses a variety of filtering methods depending on the data source. For example,
if one has matching DNA and RNA from a population, those called mutations present in both
bulk populations could preferentially be used. Or, if one has known clonal labels (e.g. barcodes 
or T-cell receptors), a semi-supervised mode of filtering may be preferred. Examples of modes of
variant filtering are present in the `exampleVariantFiltering` folder.

### Genotyping samples

#### Step 1

First, perform a pileup on each individual bam file.

```
cd exampleProcessing

# Perform a pileup for all bam samples
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-90.st.bam processed_data/BM0820_HSC_90 16571 0 BM0820_HSC_90 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-91.st.bam processed_data/BM0820_HSC_91 16571 0 BM0820_HSC_91 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-92.st.bam processed_data/BM0820_HSC_92 16571 0 BM0820_HSC_92 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-95.st.bam processed_data/BM0820_HSC_95 16571 0 BM0820_HSC_95 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-96.st.bam processed_data/BM0820_HSC_96 16571 0 BM0820_HSC_96 0
# one can use a bash loop to make this more elegant depending on file names
```

Here, the parameters descriptions in the order that they are supplied to the `01_pileup_counts.py` script. 

```
1 - Bam file for genotyping. This file should contain only mtDNA reads (extract them with samtools view)
2 - Out prefix. Filepath / basename for the single-sample
3 - max BP. Maximum length of the mtDNA genome. For hg19, 16571. For rCRS, 16569.
4 - minimum base quality. Nucleotides will only be counted if they exceed this base quality value.
5 - Sample. unique cell/sample ID for the current sample. Will be used later on when data are merged. 
6 - minimumalignment quality. Nucleotides on reads will only be counted if the read's alignment quality exceeds this value. 
```


#### Step 2

Then, merge the raw counts together:

```
sh 02_merge_pileup_counts.sh processed_data BM0820_example
```

Here, the first argument is the directory of the single samples and the second argument will be the basename of the combined
sample counts. 

#### Step 3

Finally, import the data into a [MultiAssayExperiment](http://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).

```
Rscript 03_makeRDS.R processed_data
```

where the first and only parameter is the directory that has the processed, combined data file from step 2. 

The MultiAssayExperiment is the basic unit of processed data for all of our downstream analyses. One can easily
read in the RDS file from this script and perform all downstream analyses.

For different ways of filtering variants from a MAE, see the `exampleVariantFiltering` folder. 

For an example downstream heatmap generation, see the `exampleDownstream` folder.

### Note

Code is distributed as-is noting that these scripts were used for real analyses in the manuscript. 
An optimized, efficient, transferrable pipeline is currently under development and will be available
in mid-2019.

<br><br>