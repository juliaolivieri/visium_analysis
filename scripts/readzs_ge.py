import argparse 
import pandas as pd

def get_args():
  parser = argparse.ArgumentParser(description='parse readzs counts')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--thresh', type=int, help="number of counts per window to threshold by")

  args = parser.parse_args()
  return args


def main():
  args = get_args()

  # read in dataset paths
  samples = pd.read_csv("../../notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  row = samples.loc[args.dataname]

  outpath = "../output/readzs_ge/"

  # load counts file
  df = pd.read_csv(row["readzs_counts"],sep="\t",header=None)
  df.rename(columns={0 : "cell_id",1  : "chr", 2 : "pos", 3 : "strand", 4 : "count", 5 : "sample", 6 : "window"},inplace=True)
  print("1 df",df.shape)

  # read in ReadZS zscore file
  readzs = pd.read_csv(row["readzs_vals"],sep="\t")
  print("readzs",readzs.shape)

  # only include counts that a z score was calculated for
  df = df[df["window"].isin(set(readzs["window"].unique()))]
  print("2 df",df.shape)

  # get the number of counts of each window per spot
  df["id"] = df["window"] + "_" + df["cell_id"]
  df["window_count"] = df["id"].map(df.groupby("id")["count"].sum())
  df.drop("count",axis=1,inplace=True)
  
  # now that we've summed values, we can drop duplicate cell/gene combos
  df = df.drop_duplicates(["id"])

  # merge with metadata
  meta = pd.read_csv(row["metadata"],sep="\t")
  df = df.merge(meta,on="cell_id")

  # subset to windows with number nonzero spots > thresh
  vc = df["window"].value_counts()
  df = df[df["window"].isin(vc[vc > args.thresh].index)]
  print("3 df",df.shape)

  # find the count across all windows for each cell
  df["cell_count"] = df["cell_id"].map(df.groupby("cell_id")["window_count"].sum())

  # find the fraction of each window in each cell
  df["frac_count"] = df["window_count"]/df["cell_count"]

  
  # save files
  df.to_csv("{}{}_readzs_ge_{}.tsv".format(outpath,args.dataname, args.thresh),sep="\t",index=False)
  df.to_parquet("{}{}_readzs_ge_{}.pq".format(outpath,args.dataname, args.thresh))

main()
