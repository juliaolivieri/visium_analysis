import argparse
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def get_args():
  parser = argparse.ArgumentParser(description='make spatial plots')
  parser.add_argument('--dataname', help="name of dataset")
  parser.add_argument('--score', help="first score to plot",choices=["ge","SpliZ","ReadZS","ReadZS_ge", "ReadZS_resid","ReadZS_norm","SpliZ_resid","SpliZ_norm","ReadZS_ge_norm","ge_norm"])
  parser.add_argument('--score2', help="second score to plot",choices=["ge","SpliZ","ReadZS","ReadZS_ge", "ReadZS_resid","ReadZS_norm","SpliZ_resid","SpliZ_norm","ReadZS_ge_norm","ge_norm"])
  parser.add_argument('--window_file', help="file with no header, window on each line, saying which genes/windows to plot")
  args = parser.parse_args()
  return args

def main():
  args = get_args()
  print("got args")
  outpath = "../output/plot_val_ge/"

  windows = list(pd.read_csv(args.window_file,header=None)[0]) 
  print(windows)

  scores = pd.read_csv("../../notebooks/output/make_samplesheet/scores.csv",index_col=0)
  srow = scores.loc[args.score]
  srow2 = scores.loc[args.score2]
  
  samples = pd.read_csv("../../notebooks/output/make_samplesheet/spatial.csv",index_col = 0)

  # samples
  row = samples.loc[args.dataname]
  
  
  # read in dataframes
  if row[srow["valname"]].endswith(".pq"):
    spliz_df = pd.read_parquet(row[srow["valname"]])
  else:
    spliz_df = pd.read_csv(row[srow["valname"]],sep="\t")
    
  
    
  if row[srow2["valname"]].endswith(".pq"):
    spliz_df2 = pd.read_parquet(row[srow2["valname"]])
  else:
    spliz_df2 = pd.read_csv(row[srow2["valname"]],sep="\t")
  print("read dfs")
  
  spliz_df["gene_cell"] = spliz_df[srow["genecol"]] + "_" + spliz_df[srow["cellid"]]
  spliz_df2["gene_cell"] = spliz_df2[srow2["genecol"]] + "_" + spliz_df2[srow2["cellid"]]

  thresh = 5
  xcol = "plot_xcoord"
  ycol = "plot_ycoord"
  alpha = 0.2

  # load in image
  Image.MAX_IMAGE_PIXELS = 693068558
  im = Image.open(row["image"])
  graydf = spliz_df.drop_duplicates(srow["cellid"])


  thresh = 0
  upper = 1
  lower = 0
  palette0 = "viridis"
  palette1 = "viridis"
  legval = True
  
  for gene in windows:
    print(gene)
    gene_df = spliz_df[spliz_df[srow["genecol"]] == gene]
    gene_df2 = spliz_df2[spliz_df2[srow2["genecol"]] == gene]
  
    if gene_df.shape[0] > thresh:
      
      gene_df[srow["col"] + "_quant"] = 1
  
      gene_df.loc[gene_df[srow["col"]] > gene_df[srow["col"]].quantile(.25),srow["col"] + "_quant"] = 2
      gene_df.loc[gene_df[srow["col"]] > gene_df[srow["col"]].quantile(.5),srow["col"] + "_quant"] = 3
      gene_df.loc[gene_df[srow["col"]] > gene_df[srow["col"]].quantile(.75),srow["col"] + "_quant"] = 4
      
      gene_df2[srow2["col"] + "_quant"] = 1
  
      gene_df2.loc[gene_df2[srow2["col"]] > gene_df2[srow2["col"]].quantile(.25),srow2["col"] + "_quant"] = 2
      gene_df2.loc[gene_df2[srow2["col"]] > gene_df2[srow2["col"]].quantile(.5),srow2["col"] + "_quant"] = 3
      gene_df2.loc[gene_df2[srow2["col"]] > gene_df2[srow2["col"]].quantile(.75),srow2["col"] + "_quant"] = 4
      
      gene_df[srow2["col"]] = gene_df["gene_cell"].map({x : y for x, y in zip(gene_df2["gene_cell"],gene_df2[srow2["col"]])})
      gene_df[srow2["col"] + "_quant"] = gene_df["gene_cell"].map({x : y for x, y in zip(gene_df2["gene_cell"],gene_df2[srow2["col"] + "_quant"])})
      gene_df["diff"] = gene_df[srow["col"]] - gene_df[srow2["col"]]
      gene_df["diff_quant"] = gene_df[srow["col"] + "_quant"] - gene_df[srow2["col"] + "_quant"]
      gene_df["diff2"] = (5 - gene_df[srow["col"]]) - gene_df[srow2["col"]]
      gene_df["diff2_quant"] = (5 - gene_df[srow["col"] + "_quant"]) - gene_df[srow2["col"] + "_quant"]
      for suff in ["","_quant"]:
        fig, axs = plt.subplots(1,3, figsize=(15,5))
  
        sns.scatterplot(ax=axs[1],data = graydf, x=xcol, y = ycol,color="gray",alpha = alpha,size=3,linewidth=0,legend=False)
        sns.scatterplot(ax=axs[1],data = gene_df2, x = xcol, y = ycol, hue = srow2["col"] + suff,size=3,linewidth=0,palette=palette0,legend=False)
        axs[1].set(xlabel=None)
        axs[1].set(ylabel=None)      
        axs[1].axes.xaxis.set_visible(False)
        axs[1].axes.yaxis.set_visible(False)
        sns.scatterplot(ax=axs[2],data = graydf, x=xcol, y = ycol,color="gray",alpha = alpha,size=3,linewidth=0,legend=False)
        sns.scatterplot(ax=axs[2],data = gene_df, x = xcol, y = ycol, hue = srow["col"] + suff,size=3,linewidth=0,palette=palette1,legend=False)    
        axs[2].set(xlabel=None)
        axs[2].set(ylabel=None) 
        
        
        axs[2].axes.xaxis.set_visible(False)
        axs[2].axes.yaxis.set_visible(False)
        axs[0].axes.xaxis.set_visible(False)
        axs[0].axes.yaxis.set_visible(False)
        axs[1].set_title(args.score2 + "\n")
        axs[2].set_title(args.score + "\n")
        axs[2].set_xlim([min(axs[1].get_xlim()[0],axs[2].get_xlim()[0]),max(axs[1].get_xlim()[1],axs[2].get_xlim()[1])])
        axs[2].set_ylim([min(axs[1].get_ylim()[0],axs[2].get_ylim()[0]),max(axs[1].get_ylim()[1],axs[2].get_ylim()[1])])
        axs[1].set_xlim([min(axs[1].get_xlim()[0],axs[2].get_xlim()[0]),max(axs[1].get_xlim()[1],axs[2].get_xlim()[1])])
        axs[1].set_ylim([min(axs[1].get_ylim()[0],axs[2].get_ylim()[0]),max(axs[1].get_ylim()[1],axs[2].get_ylim()[1])])
        axs[0].imshow(im)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
        try:
          axs[0].set_title("{} {} {} \n{}\n".format(args.dataname,suff,gene,name_dict[gene]))
        except:
          axs[0].set_title("{} {} {}\n".format(args.dataname,suff,gene))
  
        plt.savefig("{}{}_{}_{}_{}{}.png".format(outpath,args.dataname,gene,args.score,args.score2,suff),bbox_inches="tight")
        plt.show()
main()
