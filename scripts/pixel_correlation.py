import argparse
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

def get_args():
    parser = argparse.ArgumentParser(description='correlate pixel value with score')
    parser.add_argument('--dataname', help='dataset to correlate')
    parser.add_argument('--score',help='score to use')
    parser.add_argument('--outpath', help='path to save output',default="/scratch/groups/horence/JuliaO/single_cell/spatial_kmers/scripts/output/pixel_correlation/")
    parser.add_argument('--thresh', type=int, help='number of non-nans to require per gene')

    args = parser.parse_args()
    return args

def main():

  # get command line arguments
  args = get_args()

  # load in samples and scores
  samples = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]
  
  scores = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/scores.csv",index_col=0)
  srow = scores.loc[args.score]

  # load in dataframe
  spliz_df = pd.read_csv(row[srow["valname"]],sep="\t") 

  # subsetting to spots where pixval isn't NA
  spliz_df = spliz_df[~spliz_df["pixval"].isna()]

  vc = spliz_df[srow["genecol"]].value_counts()

  corr_cols = ["pixval", "pixquant"]
  out = {srow["genecol"] : [], "num_spots" : []}
  for corr in corr_cols:
    out["corr_" + corr] = []
    out["pval_" + corr] = []
  count = 0

  # looping over all the genes that have enough spots
  for gene in tqdm(vc[vc > args.thresh].index):
    count += 1

    # subsetting to that gene in the dataframe
    gene_df = spliz_df[spliz_df[srow["genecol"]] == gene]
  
    out[srow["genecol"]].append(gene)
    for corr in corr_cols:

      # perform the correlation
      corr_out = spearmanr(gene_df[srow["col"]],gene_df[corr])

      # save output
      out["corr_" + corr].append(corr_out.correlation)
      out["pval_" + corr].append(corr_out.pvalue)
  
    # number of spots in that gene
    out["num_spots"].append(gene_df.shape[0])
  
  out = pd.DataFrame.from_dict(out)
  out = out.dropna()
  for corr in corr_cols:
    try: 
      out["pval_{}_adj".format(corr)] = multipletests(out["pval_{}".format(corr)],method="fdr_bh")[1]
    except Exception as e:
      print(e)
      out["pval_{}_adj".format(corr)] = 1
  out.sort_values("pval_{}_adj".format(corr_cols[0]))
  out.sort_values("pval_{}_adj".format(corr_cols[0])).to_csv("{}{}_{}_{}.tsv".format(args.outpath,args.dataname,srow["col"],args.thresh),sep="\t",index=False)


main()
