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

Example call to spaceranger count with input descriptions from https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count#count

```
$ cd /home/jdoe/runs
$ spaceranger count --id=sample345 \ #Output directory
                   --transcriptome=/opt/refdata/GRCh38-2020-A \ #Path to Reference
                   --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ #Path to FASTQs
                   --sample=mysample \ #Sample name from FASTQ filename
                   --image=/home/jdoe/runs/images/sample345.tiff \ #Path to brightfield image 
                   --slide=V19J01-123 \ #Slide ID
                   --area=A1 \ #Capture area
                   --localcores=8 \ #Allowed cores in localmode
                   --localmem=64 #Allowed memory (GB) in localmode 
```

Example call to spaceranger in this repo: [`run_spaceranger.sh`](scripts/submission_scripts/run_spaceranger.sh)

The Reference was downloaded using this command:

```
curl -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

The `slide` and `area` values should be available from wherever you downloaded the data. For example, see "visium gene expression slideserial number" in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5420749.

Once you run spaceranger, the required files from [Step 1](#Step-1:-Download-data) should be available in `outs/filtered_feature_bc_matrix/` and `outs/spatial/`.

The 10X team are very helpful with debugging, so if you run into errors definitely email them.

## Step 3: Transform metadata into form required for SpliZ/ReadZS

Use the notebook [`visium_meta_clean.ipynb`](notebooks/visium_meta_clean.ipynb) to create the metadata file. Change the following at the top of the notebook to match your data:
* `dataname`: Name to use when saving the output file
* `cer_loc_name`: Path to the `tissue_positions_list.csv` file discussed in Step 1.
* `im_path`: Path to the tif image
* `filt_bc_name`: path to the filtered gene expression matrix discussed in Step 1.

Then run all the cells.

Note: sometimes I've found the images to be rotated according to what I'd expect. For example, you can see that [this pixel map](notebooks/output/visium_meta/p20190-s003_3_BrainMetastasis_pixval_mismatch.png) of a brain metastasis slide seems to be mismatched - the image is flipped relative to what we're plotting, so we're extracting the wrong pixel values from the image. You can flip around the blur matrix and transpose it until it [matches up with the histology](notebooks/output/visium_meta/p20190-s003_3_BrainMetastasis_pixval.png). You can also change the `plot_xcoord` and `plot_ycoord` columns to flip the image around so it matches the histology. I haven't found a programmatic way of knowing when this flip is happening.

Outputs:
* [`notebooks/output/visium_meta/<dataname>_blur.png`](notebooks/output/visium_meta/V1_Mouse_Brain_Sagittal_Posterior_blur.png): Plot of the histology image blurred to the specified level (default: 70)
* [`notebooks/output/visium_meta/<dataname>_pixquant.png`](notebooks/output/visium_meta/V1_Mouse_Brain_Sagittal_Posterior_pixquant.png): Plot of the pixel values from the blurred image, quantiled
* [`notebooks/output/visium_meta/<dataname>_pixval.png`](notebooks/output/visium_meta/V1_Mouse_Brain_Sagittal_Posterior_pixval.png): Plot of the raw pixel values from the blurred image
* [`notebooks/output/visium_meta/meta_<dataname>.tsv`](notebooks/output/visium_meta/meta_V1_Mouse_Brain_Sagittal_Posterior.tsv): The metadata tsv file

Columns in metadata file (some columns from `tissue_positions_list.csv` directly, see documentation here: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images):
* `barcode`: Barcode that can be mapped to BAM barcode
* `in_tissue`: 1 if in the tissue, 0 if out. All should be 1 because we're subsetting to only spots in the tissue here.
* `array_row`: The row coordinate of the spot in the array from 0 to 77. The array has 78 rows.
* `array_col`: The column coordinate of the spot in the array. In order to express the orange crate arrangement of the spots, this column index uses even numbers from 0 to 126 for even rows, and odd numbers from 1 to 127 for odd rows. Notice then that each row (even or odd) has 64 spots.
* `xcoord`: The row pixel coordinate of the center of the spot in the full resolution image.
* `ycoord`: The column pixel coordinate of the center of the spot in the full resolution image.
* `cell_id`: `<dataname>_<barcode[:-2]` (cutting off the `-1` at the end of the barcode; this should match up with names of cells in ReadZS/SpliZ)
* `plot_xcoord`: `ycoord` (use for plotting to align with histology image)
* `plot_ycoord`: `-xcoord` (use for plotting to align with histology image)
* `pixval`: pixel value at the given spot (with blur)
* `pixquant`: pixel quantile at the given spot (with blur; default is 10 quantiles)

## Step 4: Run SpliZ/ReadZS

### ReadZS

Visium data can be run using the main branch of the [ReadZS pipeline](https://github.com/salzmanlab/ReadZS) without modification. An example [config file](nextflow_inputs/visium_readzs.config), [samplesheet](nextflow_inputs/samplesheet_readzs.csv), and [bash script](nextflow_inputs/run_readzs.sh) are provided. 

In the bash script, change `runName` to whatever you want. In the samplesheet, your sample name and its BAM (I haven't run with multiple visium samples with ReadZS, but I suppose that you could). In the config file, include the path to the samplesheet, metadatafile described in Step 3, and other [input files required for ReadZS](https://github.com/salzmanlab/ReadZS#input-arguments). 

### SpliZ

## Step 5: Extract gene expression values for comparison

### By gene

### By window

## Step 6: Normalize data

## Step 7: Identify spatial patterns

### Pixel correlation

### Ising metric

## Step 8: Plot genes of interest

## References

ReadZS

SpliZ
