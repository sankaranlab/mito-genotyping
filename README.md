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

### Execution

First, perform a pileup on each individual bam file.

```
cd exampleProcessing

# Perform a pileup for all bam samples
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-90.st.bam processed_data/BM0820_HSC_90 chrM 16571 0 BM0820_HSC_90 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-91.st.bam processed_data/BM0820_HSC_91 chrM 16571 0 BM0820_HSC_91 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-92.st.bam processed_data/BM0820_HSC_92 chrM 16571 0 BM0820_HSC_92 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-95.st.bam processed_data/BM0820_HSC_95 chrM 16571 0 BM0820_HSC_95 0
python 01_pileup_counts.py raw_data/mt_singles-BM0828-HSC-frozen-151027-96.st.bam processed_data/BM0820_HSC_96 chrM 16571 0 BM0820_HSC_96 0
# one can use a bash loop to make this more elegant depending on file names
```

Then, merge the raw counts together:

```
sh 02_merge_pileup_counts.sh processed_data BM0820_example
```

Next, import the data into a [MultiAssayExperiment](http://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html).

```
Rscript 03_
```

The MultiAssayExperiment is the basic unit of processed data for all of our downstream analyses.

For different ways of filtering variants from a MAE, see the `exampleVariantFiltering` folder. 

For an example downstream heatmap generation, see the `exampleDownstream` folder.

### Note

Code is distributed as-is noting that these scripts were used for real analyses in the manuscript. 
An optimized, efficient, transferrable pipeline is currently under development and will be available
in mid-2019.

<br><br>