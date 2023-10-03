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

# Set cluster specifications, check args -----------------------------------
usage = "Usage: Rscript script.R  <prefix> <count_data> \n"
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Need at least one string")
} else {
  prefix <- args[1]
}


# For debug

# set the file prefix
prefix <- "ML_Macseq2D_MOA"
count_data <- "data/dat"

prefix
count_data

# Load the relevant libraries -----------------------------
library(Seurat)
library(edgeR)
library(data.table)
library(DT)
library(here)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyverse)
library(readxl)

# Functions ------------------------------------------------
create_dt <- function(inputdata){
  
  DT::datatable(inputdata, 
                extensions = c('Buttons', 'Scroller', "FixedColumns",'ColReorder',"Responsive"),
                rownames = FALSE, 
                options = list(dom = 'Blfrts',
                               columnDefs = list(list(className = 'dt-center',targets="_all")),
                               lengthMenu = list(c(5, -1), c('5', 'All')),
                               buttons = list(list(extend = 'colvis'),
                                              list(extend = 'csv',
                                                   filename = gsub(" ", "", paste(prefix,
                                                                                  "_AssayQC"))),
                                              list(extend = 'excel',
                                                   filename = gsub(" ", "", paste(prefix,
                                                                                  "_AssayQC")),
                                                   title = NULL)),
                               rownames = FALSE,
                               colReorder = TRUE,
                               searching = TRUE,
                               scrollX = TRUE))
}

# Read in data ---------------------------------------------------
# Preparing the bloody annotations
annotation <- 
  readr::read_csv(here::here("data/annotation/Updated_Macseq_3D_RNAseq.csv"),
                  col_names = TRUE) %>% 
  select(Barcode, row, column, target_coordinates.2, cmp, category, pathway) %>% 
  mutate(meta.annot = paste(cmp, category, pathway, sep = "-")) %>% 
  mutate(cmp.cat = paste0(cmp, "-", category)) 

## Compiling Count Reads for Seurat ----------------------------
#Read in RNAseq data
RNAseq_obj <- 
  Seurat::Read10X(here::here(count_data)) %>% 
  Seurat::CreateSeuratObject(counts = ., project = "ML_Macseq2D_MOA", min.cells = 3, min.features = 200)

RNAseq_obj[["percent.mt"]] <- 
  Seurat::PercentageFeatureSet(RNAseq_obj, pattern = "^MT-")

# Output QC data to csv ----------------------------------------
### To future whoever is using is, sorry I'm too lazy to cleanup the columns
QC.data <- 
  RNAseq_obj@meta.data %>% 
  mutate(., 
         Barcode = rownames(.)) %>% 
  tibble::tibble() %>% 
  left_join(.,annotation) %>% 
  tidyr::pivot_longer(., cols=c(2:4), names_to = 'meta.variable') %>% 
  dplyr::mutate(meta.variable = factor(meta.variable, levels = c("nCount_RNA",
                                                                 "nFeature_RNA",
                                                                 "percent.mt")),
                category = factor(category, levels = c("untreated",
                                                       "very low",
                                                       "low",
                                                       "med",
                                                       "high"))) %>% 
  
  mutate(meta.variable.2 = case_when(meta.variable == "nCount_RNA" ~ "Total mapped reads", 
                                     meta.variable == "nFeature_RNA" ~ "Total genes detected",
                                     TRUE ~ "Percent.Mito.genes")) 

QC.data %>% 
  filter(category == "untreated") %>% 
  write.csv(., paste0(here::here(), "/output/Macseq 2D MOA plate untreated samples sanitised for paper figure.csv"), row.names = FALSE)

create_dt(QC.data)


## Visualize QC metrics ------------------------------------------------------

###Overall RNAseq QC for the RNAseq output

QC.data %>%
  mutate(meta.variable.2 = factor(meta.variable.2, levels = c("Total mapped reads", 
                                                              "Total genes detected", 
                                                              "Percent.Mito.genes"))) %>% 
  ggplot(., 
         aes(x = category, y= value, fill = meta.variable.2)) +
  geom_violin(trim=FALSE)+
  scale_fill_manual(values = c("Total mapped reads" = "#58508d",
                               "Total genes detected" = "#ffa600",
                               "Percent.Mito.genes" = "#ff6361")) +
  geom_boxplot(width=0.1) +
  labs(title="Mac-seq 2D MOA plate QC metrics",
       x = "Treatment Dosage")+
  facet_wrap(.~meta.variable.2, scales = "free", ncol = 4)+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.position = "none",
        axis.text.x = element_text(family="Helvetica", face="bold", size = 10, hjust
                                   = 1, vjust = 1.1, angle = 45,
                                   margin=margin(10,0,0,0), colour = "black"))



ggsave(paste0(here::here(),"/output/KS meeting 2D Mac-seq pilot RNAseq QC metrics by MOA drug dosages.png"),width=10 ,height=5,units = "in", dpi=600)


### RNAseq QC for the control wells
QC.data %>% 
  filter(category == "untreated") %>% 
  mutate(meta.variable.2 = factor(meta.variable.2, levels = c("Total mapped reads", 
                                                              "Total genes detected", 
                                                              "Percent.Mito.genes"))) %>% 
  mutate(cmp = factor(cmp, levels = c("Media", "DMSO"))) %>% 
  ggplot(., 
         aes(x = cmp, y = value)) +
  geom_violin(trim = FALSE)+
  geom_boxplot(width = 0.1) +
  theme_bw()+
  facet_wrap(.~meta.variable.2, scales = "free")+
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(family="Helvetica", face="bold", size = 10, hjust
                                   = 1, vjust = 1.1, angle = 45,
                                   margin=margin(10,0,0,0), colour = "black"),
        axis.text.y = element_text(family="Helvetica", face="bold", size = 10,
                                   colour = "black"))


# ggsave("Fig 3 a 2D Mac-seq pilot QC metrics in negative control.png",width=10 ,height=5,units = "in", dpi=600)


### RNAseq read counts platmap by scaling ---------------------------------------

# Look at read counts of individual wells to scan for edge effects or extreme reads by scaling data against min/max

# make plots using Min-Max read counts
heapmap <-
  QC.data %>% 
  filter(meta.variable.2 == "Total mapped reads") %>% 
  mutate(Norm_value=round(((value-min(value))/(max(value)-min(value))), 2))

platetools::raw_map(data =heapmap$Norm_value,
                    well = heapmap$target_coordinates.2,
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


ggsave(paste0(here::here(),"/output/KS meeting 3D Mac-seq pilot RNAseq QC metrics normalised  total mapped reads.png"),width=10 ,height=5,units = "in", dpi=600)

## Data filter -----------------------------------------------------------------

# The "cells" used by Seurat packages in this MACseq instance refers to our treatment wells in the plate. They are filtered based on Mito content, feature numbers etc distribution.
# I filtered cells that have unique feature counts over 12500 or less than 5000, and cells that have > 20 % mitochondrial counts. 
# 
# To future self and whoever runs this: also bring positive ctrl violin plot to check the threshold. Things like death agents, then you have to allow higher percentage of mitochondria genes.

QC.Scatter.filtred <- 
  tibble(RNAseq_obj@meta.data) %>% 
  filter(nFeature_RNA > 5000 & nFeature_RNA < 12500 & percent.mt < 20) %>% 
  rename("Total mapped reads" = 1, 
         "Total genes detected" = 2,
         "Percent.Mito.genes" = 3) %>% 
  GGally::ggpairs(., 
                  columns = 2:4, 
                  aes(alpha = 0.5),
                  lower = list(continuous = "smooth")) +
  labs(title="Mac-seq pilot QC scatter plot")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

QC.Scatter.filtred
