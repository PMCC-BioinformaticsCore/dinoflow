# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# edgeR_script.R
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Authors: Mark Li, D. Ann Onuselogu
# Date of creation: 4 October 2023

## PROCESS: DIFF_GENE_EXPRESSION

# setwd("/home/aonuselogu/")

# INPUT

usage = "Usage: Rscript script.R  <seurat_obj_file_path> <annotation_file_path> <BCVplot_file_path> <output_csv_dir> <output_plot_smear_dir> \n"
args = commandArgs(trailingOnly = TRUE)

# For debugging
# seurat_obj_file_path <- "./seurat_obj.RDS"
# annotation_file_path <- "./PMC154_MACseq_MOA_annotation_modified.csv"
# BCVplot_file_path <- "./QC_out/BCV_plot.png"
# output_csv_dir <- "./CSV_out"
# output_plot_smear_dir <- "./smearPlots_out"


if (length(args) != 5) {
  stop("Oop! Incorrect number of input args. Please provide <seurat_obj_file_path> <annotation_file_path> <BCVplot_file_path> <output_csv_dir> <output_plot_smear_dir>")
} else {
  seurat_obj_file_path <- args[1]
  annotation_file_path <- args[2]
  BCVplot_file_path <- args[3]
  output_csv_dir <- args[4]
  output_plot_smear_dir <- args[5]
}


# OUTPUT

# BCV_plot.png < plotBCV output
# *_table.csv < table of genes
# _smearPlots.png < plotSmear output for each KS iteration (drug/dosage?)


# SCRIPT

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load Libraries

library(Seurat)
library(edgeR)
library(tidyverse)

# Make directories 
system(paste0("mkdir ",output_csv_dir))
system(paste0("mkdir ",output_plot_smear_dir))

# read in the object

RNAseq_obj <- readRDS(seurat_obj_file_path)
annotation <- read.csv(annotation_file_path)

# reshaping the dataset to fit edgeR
data <- GetAssayData(object = RNAseq_obj, slot = "counts")

# sort it so the order of data is checked
data <- as.data.frame(data) %>% select(sort(names(.)))

# filter data not sequenced

dataGroups <- 
  annotation %>%
  select(Barcode, cmp.cat) %>% 
  #make sure you filter the grouping factors 
  #since not all wells will be successful seq
  filter(Barcode %in% colnames(data)) %>% 
  arrange(Barcode) %>% 
  select(cmp.cat)

# Create the edgeR object

dt <- DGEList(counts = data, group = factor(dataGroups$cmp.cat))

## originally printed to console
# dt$samples

## Filter dataset
## dataset reduced from 26940 to 1746

## originally printed to console
# dim(dt)
# head(cpm(dt))
# apply(dt$counts, 2, sum) 

keep <- rowSums(cpm(dt)>100) >= 2
dt <- dt[keep,]

dt$samples$lib.size <- colSums(dt$counts)
dt$samples

dt <- calcNormFactors(dt)

## Plotting D1 --------------------------
# plotMDS(vis.dt, method="bcv", col=as.numeric(vis.dt$samples$group))
# legend("bottomleft", as.character(unique(vis.dt$samples$group)), col=1:3, pch=20)

d1 <- estimateCommonDisp(dt, verbose=T)
# names(d1)

d1 <- estimateTagwiseDisp(d1)
# names(d1)

# A generalized linear model (glm) fit to estimate dispersion
# 
# design.mat <- model.matrix(~ 0 + dt$samples$group)
# colnames(design.mat) <- levels(dt$samples$group)
# d2 <- estimateGLMCommonDisp(dt,design.mat)
# d2 <- estimateGLMTrendedDisp(d2,design.mat, method="bin.spline")
# # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
# d2 <- estimateGLMTagwiseDisp(d2,design.mat)
# plotBCV(d2)

## Dealing with D2 ----------------------

tag <- list()
sum <- list()

#KS <- c(hits$cmp.cat, "OLIGOMYCIN A-high", "Tunicamycin-high")
KS <- unique(dataGroups$cmp.cat)

# Outputs -----------------------------

## saving plot d1 - BCV plot
png(filename=BCVplot_file_path, height=700, width=800, units="px")
plotBCV(d1)
dev.off()

## saving the iterative viz plots
for (i in seq_along(KS)){
  
  et12 <- exactTest(d1, pair=c("DMSO-untreated",KS[i]))
  
  tag[[KS[i]]] <- topTags(et12, n=500)$table %>% 
    rownames_to_column(.,"genes")
  
  de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
  
  sum[[KS[i]]] <- tibble(summary(de1))%>% 
    rownames_to_column(.,"genes")
  # differentially expressed tags from the naive method in d1
  de1tags12 <- rownames(d1)[as.logical(de1)] 
  
  ## OUTPUT VIZ (iterative)
  png(filename=paste0(output_plot_smear_dir,"/",KS[i],"_smearPlots.png"), height=700, width=800, units="px")
  plotSmear(et12, de.tags=de1tags12)
  abline(h = c(-2, 2), col = "blue")
  dev.off()
  
}

## saving top 500 DEGs
tag <- data.table::rbindlist(tag,idcol=TRUE) 
write.csv(tag, paste0(output_csv_dir,"/top 500 differentially expressed genes for all treatments_table.csv"),
          row.names=F)

## saving summary of DEGs for all treatments
sum <- data.table::rbindlist(sum,idcol=TRUE) 
write.csv(sum, paste0(output_csv_dir,"/summary differentially expressed genes for all treatments_table.csv"),
          row.names=F)
