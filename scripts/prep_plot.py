import argparse
from os.path import exists
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description='save residual values as separate score')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--window_file', help="file with no header, window on each line, saying which windows to plot")

  args = parser.parse_args()
  return args

def create_zdf(row, srow, temp_path):
  # read in z scores for all
  zdf = pd.read_csv(row[srow["valname"]],sep="\t")

  # randomize and then rank so rank ties are broken evenly
  zdf = zdf.sample(frac=1)
  zdf["rank"] = zdf.groupby("window")["z_scaled"].rank("first", ascending=False)

  # assign quantile of score per window
  qs = [0.25,0.5,0.75]
  for q in qs:
    print(q)
    zdf["{}_quant".format(q)] = zdf["window"].map(zdf.groupby("window")["rank"].quantile(q))

  zdf["quant"] = 1
  for i in range(len(qs)):
    zdf.loc[zdf["rank"] > zdf["{}_quant".format(qs[i])],"quant"] = i + 2
  zdf["window_quant"] = zdf["window"] + "_" + zdf["quant"].astype(str)
  print("window_quant")
  zdf["sum_counts_per_window_per_ont"] = zdf["window_quant"].map(zdf.groupby("window_quant")["count"].sum())
  print("sum_counts_per_window_per_ont")
  zdf["med_counts_per_window_per_ont"] = zdf["window_quant"].map(zdf.groupby("window_quant")["count"].median())
  print("med_counts_per_window_per_ont")
  zdf["median_z_scaled"] = zdf["window_quant"].map(zdf.groupby("window_quant")["z_scaled"].median())
  print("median_z_scaled")
  zdf["window_cell_id"] = zdf["window"] + "_" + zdf["cell_id"]
  zdf.to_parquet(temp_path)
  return zdf

def create_zscore(row,dataname,chrom,outpath,chrompath, zdf):
  zscore = pd.read_csv("{}/zscore/{}_{}.zscore".format(row["readzs_stem"],dataname,chrom),sep="\t")
  zscore["window_cell_id"] = zscore["window"] + "_" + zscore["cell_id"]
  zscore["quant"] = zscore["window_cell_id"].map({x : y for x, y in zip(zdf["window_cell_id"],zdf["quant"])})
  zscore.to_csv(chrompath,sep="\t",index=False)
  return zscore

def make_dropped_df(zdf, ann, dropped_path):
  temp = zdf.drop_duplicates("window_quant")
  print("dropped duplicates")

  # fill in all the columns we need for plotterfile
  temp["pval"] = 0
  temp["significant"] = True
  temp["max_med"] = temp["window_quant"].map(temp.groupby("window_quant")["median_z_scaled"].max())
  temp["min_med"] = temp["window_quant"].map(temp.groupby("window_quant")["median_z_scaled"].min())
  temp["medians_range"] = temp["max_med"] - temp["min_med"]

  # add all gene names
  temp["gene"] = temp["window"].map({k : v for k, v in zip(ann["window"],ann["gene"])})
  temp["ontology"] = temp["quant"]
  temp.to_parquet(dropped_path)
  return temp


def main():
  args = get_args()
  samples = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]
  outpath = "/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/prep_plot/"
   
  score = "ReadZS"
  scores = pd.read_csv("/oak/stanford/groups/horence/JuliaO/visium_analysis/notebooks/output/make_samplesheet/scores.csv",index_col=0)
  srow = scores.loc[score]

  # load in annotated windows
  ann = pd.read_csv("{}/annotated_files/annotated_windows.file".format(row["readzs_stem"]),sep="\t")

  temp_path = "{}{}_temp.pq".format(outpath, args.dataname)

  if exists(temp_path):
    print("{} exists".format(temp_path))
    zdf = pd.read_parquet(temp_path)
  else:
    print("{} doesn't exist".format(temp_path))
    zdf = create_zdf(row, srow, temp_path)

  dropped_path = "{}{}_dropped.pq".format(outpath, args.dataname)
  if exists(dropped_path):
    print("{} exists".format(dropped_path))
    temp = pd.read_parquet(dropped_path)
  else:
    print("{} doesn't exist".format(dropped_path))
    temp = make_dropped_df(zdf, ann, dropped_path)


  windows = list(pd.read_csv(args.window_file,header=None)[0]) 
  print(windows)

  # cols for plotterfile script
  cols = ['window', 'ontology', 'sum_counts_per_window_per_ont',
             'med_counts_per_window_per_ont', 'median_z_scaled', 'pval',
                    'significant', 'medians_range', 'max_med', 'min_med', 'gene']

  for window in windows:
    print(window)
    chrom = window.split("_")[0]
    chrompath = "{}zscore/{}_{}.zscore".format(outpath,args.dataname,chrom)
    if exists(chrompath):
      print("{} exists".format(chrompath))

      zscore = pd.read_csv(chrompath,sep="\t")
    else:
      print("{} doesn't exist".format(chrompath))

      zscore = create_zscore(row,args.dataname,chrom,outpath,chrompath, zdf)
    sub = temp[temp["window"] == window][cols]
    sub.sort_values("sum_counts_per_window_per_ont").to_csv("{}pre_plotter/{}_{}.tsv".format(outpath,args.dataname,window),sep="\t", index=False)
main()
