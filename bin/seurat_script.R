# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# seurat_script.R
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Authors: Mark Li, Laura Twomey
# Date of creation: October 2023
# Version: 1.0
# Description: The aim of this analysis sequence is to prototyping the workflow for MACseq.  
# In this particular cases, 3D MCF7 cells treated with MOA plate for 24hr recovered by Macseq workflow (Cell Recovery Solution), and a parallel plate analysed via CTG assay.  <br>
# To select for treatment wells for investigation, CTG data is used to help me selecting a few candidates.  <br>
# I would like to perform RNAseq QC, and differential expression analysis of treatment wells compared to negative control DMSO. <br>

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Load the relevant libraries -----------------------------
library(Seurat)
library(data.table)
library(DT)
library(plotly)
library(tidyverse) # includes ggplot2 and dplyr
library(readxl)
library(platetools)
library(GGally)

# Functions ------------------------------------------------

# Inputs -----------------------------------
usage = "Usage: Rscript script.R  <prefix> <count_data> <mtx_file> <barcode_file> <feature_file> <output_folder> <out_RDS_file> <out_ann_file> \n"
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Oop! Incorrect number of input args. Please provide <prefix> <count_data> <mtx_file> <barcode_file> <feature_file> <output_folder> <out_RDS_file> <out_ann_file>")
} else {
  prefix <- args[1]
  ann_data <- args[2]
  mtx_file <- args[3]
  barcode_file <- args[4]
  feature_file <- args[5]
  output_folder <- args[6]
  out_RDS_file <- args[7]
  out_ann_file <- args[8]
}

# Columns in annotation file
ann_cols <- c('Barcode', 'Well_ID', 'Row', 'Column', 'Compound_ID', 'Dose_category', 'Pathway')

# # For debug
# prefix <- "2D_PMC154_ML_Macseq2D_MOA"
# ann_data <- "/home/ltwomey/data/Sequencing/2D_monolayer/MOA/Sequencing_2D_monolayer_MOA_Rep1and2_annotation_new.csv"
# mtx_file <- '/home/ltwomey/data/Sequencing/2D_monolayer/MOA/Rep1/preprocessed/matrix.mtx.gz'
# barcode_file <- '/home/ltwomey/data/Sequencing/2D_monolayer/MOA/Rep1/preprocessed/barcodes.tsv.gz'
# feature_file <- '/home/ltwomey/data/Sequencing/2D_monolayer/MOA/Rep1/preprocessed/features.tsv.gz'
# out_RDS_file <- '/home/ltwomey/seurat_obj.RDS'
# out_ann_file <- '/home/ltwomey/PMC154_MACseq_MOA_annotation_modified.csv'

prefix
ann_data
mtx_file
barcode_file
feature_file
out_RDS_file
out_ann_file

# Read in data ---------------------------------------------------
# Preparing the bloody annotations
annotation <- readr::read_csv(ann_data, col_names = TRUE) 
# Check that all columns are there before selecting them
stopifnot(sum(!ann_cols %in% colnames(annotation)) == 0)
# Now fix up the ann df
annotation <- annotation %>% 
  select(all_of(ann_cols)) %>% 
  mutate(meta.annot = paste(Compound_ID, Dose_category, Pathway, sep = "-")) %>% 
  mutate(cmp.cat = paste0(Compound_ID, "-", Dose_category)) 

## Compiling Count Reads for Seurat ----------------------------
expression.matrix <- ReadMtx(mtx = mtx_file, features = feature_file, cells = barcode_file)
seu <- CreateSeuratObject(counts = expression.matrix, project = prefix, min.cells = 3, min.features = 200)
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

# QC  -------------------------------------------------------------
QC.data <- seu@meta.data %>% 
  mutate(., Barcode = rownames(.)) %>% 
  tibble::tibble() %>% 
  left_join(.,annotation) %>% 
  tidyr::pivot_longer(., cols=c(2:4), names_to = 'meta.variable') %>% 
  dplyr::mutate(meta.variable = factor(meta.variable, levels = c("nCount_RNA",
                                                                 "nFeature_RNA",
                                                                 "percent.mt")),
                Dose_category = factor(Dose_category, levels = c("untreated",
                                                            "very low",
                                                            "low",
                                                            "med",
                                                            "high"))) %>% 
  
  mutate(meta.variable.2 = case_when(meta.variable == "nCount_RNA" ~ "Total mapped reads", 
                                     meta.variable == "nFeature_RNA" ~ "Total genes detected",
                                     TRUE ~ "Percent.Mito.genes")) %>%
  mutate(meta.variable.2 = factor(meta.variable.2, levels = c("Total mapped reads", 
                                                              "Total genes detected", 
                                                              "Percent.Mito.genes")))

QC.data_DMSO <- QC.data %>% 
  filter(Dose_category == "untreated") %>% 
  select(-c(cmp.cat, meta.annot)) 

### Overall RNAseq QC for the RNAseq output
QC_barplots_compounds <- ggplot(QC.data, aes(x = Dose_category, y = value, fill = meta.variable.2)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("Total mapped reads" = "#58508d",
                               "Total genes detected" = "#ffa600",
                               "Percent.Mito.genes" = "#ff6361")) +
  geom_boxplot(width = 0.1) +
  labs(title = paste(prefix, "QC metrics"),
       x = "Treatment Dosage")+
  facet_wrap(.~ meta.variable.2, scales = "free", ncol = 4)+
  theme_bw() +
  theme(axis.title.y = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.position = "none",
        axis.text.x = element_text(family = "Helvetica", face = "bold", size = 10, 
                                   hjust = 1, vjust = 1.1, angle = 45,
                                   margin = margin(10, 0, 0, 0), colour = "black"))


### RNAseq QC for the control wells
QC_barplots_untreated <- ggplot(QC.data_DMSO, aes(x = Compound_ID, y = value)) +
  geom_violin(trim = FALSE)+
  geom_boxplot(width = 0.1) +
  theme_bw() +
  facet_wrap(.~ meta.variable.2, scales = "free")+
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(family="Helvetica", face = "bold", size = 10, 
                                   hjust = 1, vjust = 1.1, angle = 45,
                                   margin = margin(10, 0, 0, 0), colour = "black"),
        axis.text.y = element_text(family = "Helvetica", face = "bold", size = 10,
                                   colour = "black"))

# RNAseq read counts platmap by scaling ---------------------------------------
# Look at read counts of individual wells to scan for edge effects or extreme reads by scaling data against min/max
# make plots using Min-Max read counts
QC.data %>% filter(meta.variable.2 == "Total mapped reads") %>% 
  mutate(Norm_value = round(((value-min(value))/(max(value)-min(value))), 2))

heapmap_plot <- platetools::raw_map(data = heapmap$Norm_value,
                    well = heapmap$Well_ID,
                    plate = 384) +
  ggtitle("2D Mac-seq: Normalised total mapped reads")+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = mean(heapmap$Norm_value),
                       limits = c(floor(0), 
                                  ceiling(max(heapmap$Norm_value))),
                       name = "RNA read counts") +
  theme(plot.title = element_text(size = 22, hjust = 0.05, vjust = -0.1, margin = margin(b = -5)),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4))


# Data filter -----------------------------------------------------------------
# THIS NEEDS FIXING!!!
# The "cells" used by Seurat packages in this MACseq instance refers to our treatment wells in the plate. They are filtered based on Mito content, feature numbers etc distribution.
# I filtered cells that have unique feature counts over 12500 or less than 5000, and cells that have > 20 % mitochondrial counts. 
# To future self and whoever runs this: also bring positive ctrl violin plot to check the threshold. Things like death agents, then you have to allow higher percentage of mitochondria genes.

QC.Scatter.filtered <- tibble(seu@meta.data) %>% 
  filter(nFeature_RNA > 5000 & nFeature_RNA < 12500 & percent.mt < 20) %>% 
  rename("Total mapped reads" = 1, 
         "Total genes detected" = 2,
         "Percent.Mito.genes" = 3) %>% 
  GGally::ggpairs(., 
                  columns = 2:4, 
                  aes(alpha = 0.5),
                  lower = list(continuous = "smooth")) +
  labs(title = "Mac-seq pilot QC scatter plot")+
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

# Outputs -----------------------------------------------------
write.csv(QC.data_DMSO, paste0(output_folder, "/untreated_QC.csv"), row.names = FALSE)
ggsave(QC_barplots_compounds, filename = paste0(output_folder,"/QC_pertreatmentdosage_compound_barplots.png"), width = 10, height = 5, units = "in", dpi = 600)
ggsave(QC_barplots_untreated, filename = paste0(output_folder,"/QC_pertreatmentdosage_untreated_barplots.png"), width = 10, height = 5, units = "in", dpi = 600)
gsave(heapmap_plot, filename = paste0(output_folder, "/QC_plate_heatmap.png"), width = 10 , height = 5, units = "in", dpi = 600)
gsave(QC.Scatter.filtered, filename = paste0(output_folder, "/QC_scatterplot.png"), width = 10 , height = 5, units = "in", dpi = 600)
# Save modified annotation object for edgeR
write.csv(x = annotation, file = out_ann_file)
# Save seurat object
saveRDS(seu, file = out_RDS_file)

