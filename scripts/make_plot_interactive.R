ptm <- proc.time()
#!/usr/bin/env Rscript
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if (!require("Gviz")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Gviz")
  library(Gviz)
}

library(data.table)
library(Gviz)



args <- commandArgs(TRUE)
# window <- "chr8_11013_plus"
# window <-"chr11_10983_plus"
window <- "chr17_16275_minus"
dataname <- "V1_Mouse_Brain_Sagittal_Posterior"
#dataname <- "V1_Mouse_Kidney"
input_file <- paste0("/oak/stanford/groups/horence/JuliaO/visium_analysis/scripts/output/peak_plot/",dataname,"_",window,".plotterFile")
ont_cols <- "quant"
gff_path <- "/oak/stanford/groups/horence/kaitlin/ref_files/gencode.vM26.annotation.gff3"

for (end in c("start","end")) {
  for (i in 1:2) {
    ont_cols <- unlist(strsplit(ont_cols, ", "))
    gen = "hg38"
    data <-fread(input_file, header=T)
    
    chr <- unique(data$chromosome)
    window_start <- unique(data$window_start)
    window_end <- unique(data$window_end)
    if (end == "start") {
      window_start = window_end -  ceiling((window_end - window_start)/i)  
    } else if (end == "end") {
      window_end = window_start + ceiling((window_end - window_start)/i)
    }
    
  #  if (window == "chr8_11013_plus") {
  #    window_end <- 55061000
  #  } else if (window == "chr11_10983_plus") {
  #    window_end <- 54910300
  #  } 
  
    data <- data[, ontology:=Reduce(function(...) paste(..., sep="___"), .SD[, mget(ont_cols)])]
    data <- data[, z_label := paste("Med Z = ", round(median_z_scaled, 3), sep="")]
    data <- data[, label_col := paste(ontology, z_label, metrics, sep="; ")]
    data$label_col <- gsub("; ","\n", data$label_col)
    data$label_col <- gsub("___", "\n", data$label_col)
    
    data <- data[, seqnames := chr]
    data <- data[, end := pos_round+1]
    data <- data[,c("pos_round", "end", "percent", "seqnames", "label_col")]
    names(data) <- c("start", "end", "percent", "seqnames","label_col")
    
    name_list<-unique(data$label_col)
    data$label_col <- factor(data$label_col, levels=unique(data$label_col))
    
    tracks <- split(data[, c("start","end","percent","seqnames")], data$label_col)
    plot_list = c()
    
    padder <- data.table()
    padder <- padder[, start := seq(window_start, window_end, 10)]
    padder <- padder[, end := start + 10]
    padder <- padder[, percent := 0]
    
    ind=1
    for (track in tracks){
      df <- data.frame(track)
      names(df) <-  c("start", "end","percent","seqnames")
      for (j in 1:nrow(padder)){
        target_start <- padder[j, start]
        target_end <- target_start + 10
        if (target_start %in% df$start == FALSE){
          new_row <- c(target_start, target_end, 0, chr)
          df <- rbind(df, new_row)
        }
      }
      xranges <- makeGRangesFromDataFrame(df, TRUE)
      xtrack<-DataTrack(
        xranges,
        name=name_list[ind],
        chromosome=chr,
        from=window_start,
        to=window_end,
        cex.axis=0.01, col="black", lwd=3)
      plot_list <- c(plot_list, xtrack)
      ind=ind+1
    }
    
    geneTrack <- GeneRegionTrack(
      range=gff_path,
      genome=gen,
      chromosome=chr,
      start=window_start,
      end=window_end,
      fontsize=8,panel.only=TRUE,
      options(ucscChromosomeNames=FALSE)
    )
    ranges(geneTrack)$feature <- as.character(ranges(geneTrack)$feature)
    ranges(geneTrack) <- ranges(geneTrack)[ranges(geneTrack)$feature!="transcript"]
    
    plot_list <- c(plot_list, geneTrack)
    
    outfile <- sub(".plotterFile", paste0("_peak_int_",window_start,"_",window_end,".png"), input_file)
    
    png(outfile, width = 800, height = 800)
    plotTracks(
      plot_list,
      from=window_start,
      to=window_end,
      shape="arrow",
      exonAnnotation='symbol',
      fontcolor.exon="black",
      fontcolor.title="black",
      type='a',
      legend=TRUE,
      cex.legend=.8,
      cex.title=1.5
    )
    dev.off()
    plotTracks(
      plot_list,
      from=window_start,
      to=window_end,
      shape="arrow",
      exonAnnotation='symbol',
      fontcolor.exon="black",
      fontcolor.title="black",
      type='a',
      legend=TRUE,
      cex.legend=.8,
      cex.title=1.5
    )
    
    total_time <- proc.time() - ptm
    print(total_time)
  }
}