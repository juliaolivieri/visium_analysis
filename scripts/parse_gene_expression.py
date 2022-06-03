import argparse
import pandas as pd
import scanpy as sp
from tqdm import tqdm

def get_args():
    parser = argparse.ArgumentParser(description='extract gene expression values')
    parser.add_argument('--dataname', help='dataset to extract values for')

    parser.add_argument('--thresh',type=int, help='number of spots to require for subset',default=500)
    parser.add_argument('--genes',nargs="*", help='specific genes to parse',default=[])
    parser.add_argument('--norm',action="store_true", help='normalize each cell')
    args = parser.parse_args()
    return args


def main():
  args = get_args()

  # read in dataset paths
  samples = pd.read_csv("../../notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]

  if args.norm:
    suff = ""
  else:
    suff = "_count"

  # if a list of genes are submitted, we'll only parse out those genes
  if len(args.genes) > 0:
    suff += "_" + "_".join(args.genes)
  outpath = "../output/parse_gene_expression/"

  # read in metadata
  meta_path = row["metadata"]
  meta = pd.read_csv(meta_path,sep="\t")

  # read in gene expression matrix from spaceranger
  data = sp.read_mtx(row["ge_mat"] + "matrix.mtx.gz")
  data = data.T
  features = pd.read_csv(row["ge_mat"] + "features.tsv.gz", header=None)
  barcodes = pd.read_csv(row["ge_mat"] + "barcodes.tsv.gz", header=None)
  data.var_names = features[0]
  data.obs_names = barcodes[0]
  geneex = data.to_df()

  # divide by number of counts for each cell (otherwise everything will correspond to read depth)
  if args.norm:
    geneex = geneex.div(geneex.sum(axis=1), axis=0)

  # create a table with gene expression in the given format
  out = {"gene" : [], "ensembl" : [], "frac_count" : [], "barcode" : []}

  # if we specify gene names, only keep the genes with the same names
  usecols = list(geneex.columns)
  if len(args.genes) > 0:
    usecols = [x for x in usecols if x.split("\t")[1] in args.genes]

  print("usecols",usecols)
  
  # for each gene and each cell, add a row to the output with its value
  for col in tqdm(usecols):
    for barcode in geneex.index:
  
      out["gene"].append(col.split("\t")[1])
      out["ensembl"].append(col.split("\t")[0])
      out["barcode"].append(barcode)
      out["frac_count"].append(geneex.loc[barcode,col])
  out = pd.DataFrame.from_dict(out)

  # include metadata in table
  out = out.merge(meta,on="barcode")
  out.to_csv("{}{}{}.tsv".format(outpath,args.dataname,suff),sep="\t",index=False)
  out.to_parquet("{}{}{}.pq".format(outpath,args.dataname,suff))

  # subset to only genes with greater than args.thresh nonzero values
  out["bool"] = 0
  out.loc[out["frac_count"] > 0, "bool"] = 1
  out["num_spot_nnz"] = out["gene"].map(out.groupby("gene")["bool"].sum())
  out[out["num_spot_nnz"] > args.thresh].to_parquet("{}{}_sub_{}{}.pq".format(outpath,args.dataname,args.thresh,suff))
  out[out["num_spot_nnz"] > args.thresh].to_csv("{}{}_sub_{}{}.tsv".format(outpath,args.dataname,args.thresh,suff),sep="\t",index=False)

main()
