import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description='save residual values as separate score')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--score', help="score to find residuals of",choices=["ge","SpliZ","ReadZS","ReadZS_ge", "ReadZS_resid","ReadZS_norm","SpliZ_resid","SpliZ_norm","ReadZS_ge_norm","ge_norm"])
  parser.add_argument('--thresh', type=int,help="minimum number of spots per window")
  parser.add_argument('--num_perms', type=int,help="number of permutations to perform per gene")

  args = parser.parse_args()
  return args

def main():
  outpath = "/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/ising/"
  args = get_args()


  samples = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]
  
  scores = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/scores.csv",index_col=0)
  srow = scores.loc[args.score]

  # open both dataframes
  df = pd.read_csv(row[srow["valname"]],sep="\t")

  radius = int(row["radius"])
  xcol = "plot_xcoord"
  ycol = "plot_ycoord"
  
  # create dataframe where each row is a pair of spots that are neighbors
  df["x_y"] = df[xcol].astype(str) + "_" + df[ycol].astype(str)
  temp = df.drop_duplicates("x_y")[[xcol,ycol,"x_y"]]
  out = {"x_y1": [], "x_y2" : []}

  # keep if they are neighbors
  for index, row in temp.iterrows():
    temp2 = temp[(temp[xcol] > row[xcol] - radius) & (temp[xcol] < row[xcol] + radius)& (temp[ycol] < row[ycol] + radius) & (temp[ycol] > row[ycol] - radius) & (temp["x_y"] != row["x_y"])]
    for ind2, row2 in temp2.iterrows():
      
      out["x_y1"].append(min([row["x_y"],row2["x_y"]]))
      out["x_y2"].append(max([row["x_y"],row2["x_y"]]))
  outdf = pd.DataFrame.from_dict(out)
  outdf = outdf.drop_duplicates()
  
  vc = df[srow["genecol"]].value_counts()
  
  out = {srow["genecol"] : [], "score_cont" : [], "num_pairs" : [], "perm_pval" : [], "mean_score" : []}
  for gene in tqdm(vc[vc > args.thresh].index):
    
    genedf = df[df[srow["genecol"]] == gene]
    
    # quantile values
    genedf["spin"] = -1
    genedf.loc[genedf[srow["col"]] > genedf[srow["col"]].quantile(0.5),"spin"] = 1

    spinsums2 = 0
    num_pairs = 0
    # get indicies of neighbors present for this gene
    ind1 = outdf[(outdf["x_y1"].isin(genedf["x_y"])) & (outdf["x_y2"].isin(genedf["x_y"]))]["x_y1"]
    ind2 = outdf[(outdf["x_y1"].isin(genedf["x_y"])) & (outdf["x_y2"].isin(genedf["x_y"]))]["x_y2"]
    num_pairs = ind1.shape[0]
    
    # index by these indices and then take the dot product
    dot_prod = np.dot(genedf.set_index("x_y").loc[list(ind1)][srow["col"]],genedf.set_index("x_y").loc[list(ind2)][srow["col"]])

    # permute the spot locations and take the dot product num_perm times
    perm_dot_prods = []
    for i in range(args.num_perms):
      genedf["perm"] = np.random.permutation(genedf[srow["col"]])
      perm_dot_prods.append(np.dot(genedf.set_index("x_y").loc[list(ind1)]["perm"],genedf.set_index("x_y").loc[list(ind2)]["perm"]))
    out["perm_pval"].append(len([x for x in perm_dot_prods if x > dot_prod])/args.num_perms)
    out[srow["genecol"]].append(gene)
    out["score_cont"].append(dot_prod/num_pairs)
    out["num_pairs"].append(num_pairs)
    out["mean_score"].append(genedf[srow["col"]].mean())

  out = pd.DataFrame.from_dict(out)
  out = out.sort_values("perm_pval")
  out["perm_pvals_adj"] = multipletests(out["perm_pval"],alpha=0.05,method="fdr_bh")[1]
  out.to_csv("{}{}_{}_{}_{}.tsv".format(outpath,args.dataname,args.score,args.thresh, args.num_perms),sep="\t",index=False)
  out[out["perm_pvals_adj"] < 0.05].sort_values("score_cont",ascending=False)[srow["genecol"]].to_csv("{}{}_{}_{}_{}_plot.txt".format(outpath,args.dataname,args.score,args.thresh, args.num_perms),index=False,header=None)
main()
