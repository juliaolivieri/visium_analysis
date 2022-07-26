import pandas as pd
import scanpy as sp
from tqdm import tqdm

def main():
  outpath = "../output/make_overview_table/"
  samples = pd.read_csv("../../notebooks/output/make_samplesheet/spatial.csv",index_col = 0)
  scores = pd.read_csv("../../notebooks/output/make_samplesheet/scores.csv",index_col=0)
  sub_scores = ["SpliZ", "ReadZS", "ge", "ReadZS_ge"]
  out = {"dataname" : [], "num_spots" : [], "med_reads_per_spot" : []}
  for score in sub_scores:
    out["{}_med_per_spot".format(score)] = []

  datanames = ['V1_Mouse_Brain_Sagittal_Posterior',
         'V1_Mouse_Brain_Sagittal_Posterior_Section_2',
         'V1_Mouse_Brain_Sagittal_Anterior',
         'V1_Mouse_Brain_Sagittal_Anterior_Section_2','V1_Mouse_Kidney',
          'Visium_FFPE_Human_Breast_Cancer', 'Visium_FFPE_Human_Normal_Prostate',
          'Visium_FFPE_Human_Prostate_Acinar_Cell_Carcinoma',
          'Visium_FFPE_Human_Prostate_Cancer', 'Visium_FFPE_Human_Prostate_IF',
          'Visium_FFPE_Mouse_Brain', 'Visium_FFPE_Mouse_Brain_IF',
          'Visium_FFPE_Mouse_Kidney',  'p20190_s003_3_BrainMetastasis',
          'p20190_s004_4_BrainMetastasis', 'p20218_s001_L1', 'p20218_s002_L2',
          'p20218_s003_L3', 'p20218_s004_L4']
  for dataname, row in tqdm(samples[samples["method"] == "visium"].iterrows()):
  
    if dataname in datanames:
      out["dataname"].append(dataname)
      meta = pd.read_csv(row["metadata"],sep="\t")
      out["num_spots"].append(meta["in_tissue"].sum())
      try:
        data = sp.read_mtx(row["ge_mat"] + "matrix.mtx.gz")
        out["med_reads_per_spot"].append(data.to_df().sum(axis=0).median())
      except:
        out["med_reads_per_spot"].append(np.nan)

      for score, srow in scores.iterrows():
        if score in sub_scores:
          try:
            df = pd.read_csv(row[srow["valname"]],sep="\t",usecols = [srow["cellid"],srow["genecol"]])
            out["{}_med_per_spot".format(score)].append(df.groupby(srow["cellid"])[srow["genecol"]].nunique().median())
  
          except:
            out["{}_med_per_spot".format(score)].append(np.nan)
  out = pd.DataFrame.from_dict(out)
  for k, v in out.items():
      print(k,len(v))
  out.to_csv("{}table.tsv".format(outpath),sep="\t",index=False)
main() 
