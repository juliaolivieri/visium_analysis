# Visium analysis with the SpliZ / ReadZS

## Step 1: Download data

All the Visium data I've downloaded is here: `/oak/stanford/groups/horence/JuliaO/data/visium/`. I copied the "batch download" commands from https://www.10xgenomics.com/resources/datasets/ and ran them in sherlock (separately for each sample). Key files are:

* `*.tif`: The histology image

### If BAM is available
As in https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard-1-1-0; see `/oak/stanford/groups/horence/JuliaO/data/visium/V1_Mouse_Brain_Sagittal_Posterior` for example downloaded data:
* `*possorted_genome_bam.bam`: The BAM file. 
* `spatial/tissue_positions_list.csv`: maps barcodes to spatial locations (described here: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images). 
* `filtered_feature_bc_matrix/*`: Gives gene expression values to compare against. Should contain `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`   

### If BAM is not available

As in https://www.biorxiv.org/content/10.1101/2021.08.03.455000v1.full.pdf; see `/oak/stanford/groups/horence/JuliaO/data/visium/brain_metastases` for example downloaded data:
* `*_L001_R{1,2}_001.fastq.gz`: R1 and R2 fastqs

## Step 2: Run spaceranger to get BAMs and spatial metadata

If the BAM was not available for the download, you need to run spaceranger to get the other required files. Spaceranger is downloaded on sherlock here: `/home/groups/horence/applications/spaceranger-1.3.1/`.

Example call to spaceranger in this repo: [`run_spaceranger.sh`](scripts/submissionscripts/run_spaceranger.sh)

## Step 3: Transform metadata into form required for SpliZ/ReadZS

## Step 4: Run SpliZ/ReadZS

## Step 5: Extract gene expression values for comparison

### By gene

### By window

## Step 6: Normalize data

## Step 7: Identify spatial patterns

### Pixel correlation

### Ising metric

## Step 8: Plot genes of interest

