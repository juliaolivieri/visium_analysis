import argparse
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description='perform correlations')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--method', help="score to run correlations on",choices=["ge","SpliZ","ReadZS"])
  parser.add_argument('--val_path', help="path to the file with the score")
  parser.add_argument('--num_genes', type=int,help="maximum number of genes to include")
  parser.add_argument('--thresh', type=int,help="minimum number of cells per window")

  args = parser.parse_args()
  return args

def get_mat(df,val,cell_id,genecol):
  mat = pd.DataFrame(index=df[cell_id].unique(),columns=sorted(df[genecol].unique()))
  
  for gene, gene_df in tqdm(df.groupby(genecol)):
      mat[gene] = mat.index.map({k : v for k, v in zip(gene_df[cell_id],gene_df[val])})

  return mat
      
def get_corr(mat, lower_lim):
  cols = mat.columns
  lower_lim = 20
  
  out = {"gene1" : [], "gene2" : [], "spearman_corr" : [], "spearman_pval" : [],"n" : []}
  for i in tqdm(range(len(cols) - 1)):
    for j in range(i + 1, len(cols)):
      temp = mat[[cols[i],cols[j]]].dropna()
      if (temp.shape[0] > lower_lim) & (temp[cols[i]].nunique() > 1) & (temp[cols[j]].nunique() > 1):
        x = spearmanr(temp[cols[i]],temp[cols[j]])
        out["gene1"].append(cols[i])
        out["gene2"].append(cols[j])
        out["spearman_corr"].append(x.correlation)
        out["spearman_pval"].append(x.pvalue)
        out["n"].append(temp.shape[0])
  
  outdf = pd.DataFrame.from_dict(out)
  outdf["pval_adj"] = multipletests(outdf["spearman_pval"])[1]
  outdf["pair"] = outdf["gene1"] + "_" + outdf["gene2"]
  return outdf

def main():
  args = get_args()
  samples = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]

  outpath = "/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/visium_correlations/"
  supdir = "/oak/stanford/groups/horence/JuliaO/nf-spliz-output/run_dir/visium/"

  if args.method == "ge":
    val = "frac_count"
    genecol = "gene"
    val_path = "ge_vals"
#    val_path = "/scratch/groups/horence/JuliaO/single_cell/spatial_kmers/scripts/output/parse_gene_expression/{}.pq".format(args.dataname)
    cell_id = "cell_id"
    # spliz_path = "/scratch/groups/horence/JuliaO/single_cell/spatial_kmers/scripts/output/parse_gene_expression/V1_Mouse_Kidney_count_Usp9y_Kdm5d.pq"
  
  elif args.method == "ReadZS":
    val = "z_scaled"
    genecol = "window"
    val_path = "readzs_vals"
#    val_path = "/oak/stanford/groups/horence/JuliaO/data/visium/{}/{}_all_5000.zscore".format(args.dataname,args.dataname)
    cell_id = "cell_id"
  elif args.method == "SpliZ":
    val = "scZ"
    genecol = "gene"
    val_path = "spliz_vals"
#    val_path = "{}{}/SpliZ_values/{}_sym_SVD_normdonor_S_0.1_z_0.0_b_5_r_0.01_subcol.tsv".format(supdir,args.dataname,args.dataname)
    cell_id = "cell"    

  df = pd.read_csv(row[val_path],sep="\t")
#  if args.val_path.endswith(".pq"):
#    df = pd.read_parquet(args.val_path)
#  elif args.val_path.endswith(".tsv"):
#    df = pd.read_csv(args.val_path,sep="\t")

  vc = df[genecol].value_counts()
  vc = vc.head(args.num_genes)
  df = df[df[genecol].isin(vc[vc > args.thresh].index)]

  mat = get_mat(df,val,cell_id,genecol)


  lower_lim = 20
  outdf = get_corr(mat, lower_lim)
  outdf.sort_values("spearman_pval").to_csv("{}{}_{}_{}_{}_corr.tsv".format(outpath,args.dataname, args.method, args.num_genes, args.thresh),index=False,sep="\t")



main()


