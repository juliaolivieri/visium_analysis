import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scipy.stats import norm
from scipy.stats import linregress
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description='save residual values as separate score')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--score', help="score to find residuals of",choices=["ge","SpliZ","ReadZS","ReadZS_ge"])
  parser.add_argument('--score2', help="score to regress out from previous",choices=["ge","SpliZ","ReadZS","ReadZS_ge"])
  parser.add_argument('--outname', help="path to save output")
  parser.add_argument('--thresh', type=int,help="minimum number of spots per window")

  args = parser.parse_args()
  return args

def get_res_dict(q, srow, srow2,df):
  cols = [srow["col"],srow2["col"]]
  out = {srow["genecol"] : [], "num_spots" : [], "r_squared" : []}
  
  res_map = {}
  norm_map = {}
  norm_map2 = {}

  # get residual for each window
  for window, windf in tqdm(df.groupby(srow["genecol"])):
    for col in cols:
#      qval = q
#      
#      # quantile
#      while True:
#        try:
#          windf["{}_quant".format(col)] = pd.qcut(windf[col],qval,labels=False)
#  
#          break
#        except:
#          qval -= 1

      # randomize rows so ties are broken randomly (for ranking)
      windf = windf.sample(frac=1)
      
      # rank each  value
      windf["{}_quant".format(col)] = windf[col].rank(method="first")
  
      # inverse cdf to normal
      # get uniform distribution of quantile values in (0,1) and then convert
      windf["{}_norm".format(col)] = norm.ppf((windf["{}_quant".format(col)])/(windf["{}_quant".format(col)].max() + 1))
    
    # perform regression on these values
    result = linregress(windf["{}_norm".format(srow["col"])], windf["{}_norm".format(srow2["col"])])
    out[srow["genecol"]].append(window)
    out["num_spots"].append(windf[srow["cellid"]].nunique())
    out["r_squared"].append(result.rvalue**2)

    # find prediction from regression, and residuals
    windf["predict"] = result.intercept + result.slope*windf["{}_norm".format(srow2["col"])]
    windf["res"] = windf["{}_norm".format(srow["col"])] - windf["predict"]
    norm_map = {**norm_map,**{x : y for x, y in zip(windf["gene_cell"],windf["{}_norm".format(srow["col"])])}} 
    norm_map2 = {**norm_map2,**{x : y for x, y in zip(windf["gene_cell"],windf["{}_norm".format(srow2["col"])])}} 

    res_map = {**res_map,**{x : y for x, y in zip(windf["gene_cell"],windf["res"])}}
  return res_map, norm_map, norm_map2

def main():
  args = get_args()

  # q is not used in this version
  q = 10

  # load in samples and scores
  samples = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]
  
  scores = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/scores.csv",index_col=0)
  srow = scores.loc[args.score]
  srow2 = scores.loc[args.score2]

  # open both dataframes
  print(row[srow["valname"]])
  df = pd.read_csv(row[srow["valname"]],sep="\t")
  print(row[srow2["valname"]])
  df2 = pd.read_csv(row[srow2["valname"]],sep="\t")

  # subset for number of spots (in score 1 df)
  df["num_spots"] = df[srow["genecol"]].map(df.groupby(srow["genecol"])[srow["cellid"]].nunique())
  df = df[df["num_spots"] > args.thresh]
  print("subset")

  # create column to use to map between scores
  df["gene_cell"] = df[srow["genecol"]] + "_" + df[srow["cellid"]]
  df2["gene_cell"] = df2[srow2["genecol"]] + "_" + df2[srow2["cellid"]]
  print("joined column")

  # include score2 in score1 df
  df[srow2["col"]] = df["gene_cell"].map({x : y for x, y in zip(df2["gene_cell"],df2[srow2["col"]])})
  print("include score")

  # only include rows where neither score is NA
  df = df[(~(df[srow["col"]].isna())) & (~(df[srow2["col"]].isna()))]
  print("subset to not NA")

  # get res_dict
  res_map, norm_map, norm_map2 = get_res_dict(q, srow, srow2,df)
  print("get res dict")

  df["res"] = df["gene_cell"].map(res_map)
  df["{}_norm".format(srow["col"])] = df["gene_cell"].map(norm_map)
  df["{}_norm".format(srow2["col"])] = df["gene_cell"].map(norm_map2)

  df = df[~df["res"].isna()]

  # save output
  df.to_csv(args.outname,sep="\t",index=False)

  
main()
