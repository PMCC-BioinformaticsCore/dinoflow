library(R.utils)
library(tidyverse)

#fetch the id to create samplesheets

project_name<-"nsd_pilot"
platform<-"NovaSeq"
run_id<-"230712_A01524_0152_AH5HWKDRX3"
researcher<-"Lizzy Pijpers"

metadata_dir<-paste0("/researchers/nenad.bartonicek/projects/mac-seq/metadata/",project_name)
samplesheet_dir<-paste0("/researchers/nenad.bartonicek/projects/mac-seq/sample_sheets/",project_name)

metadata_file<-paste0(metadata_dir,"/",project_name,".csv")
metadata<-read.csv(metadata_file)

#for each plate create new sample sheet
pools<-unique(metadata$Plate_ID)

dfL<-list()
for(pool in pools){
  sample_name<-paste0("Sample_",pool)
  files<-list.files(
    paste0(
      "/pipeline/Archives/",
      platform,
      "/",run_id,
      "/ProjectFolders",
      "/Project_",gsub(" ","-",researcher),
      "/",sample_name
    ),full.names=T
  )
  nLanes<-sum(grepl("R1",basename(files)))
  tempL<-list()
  for(lane in 1:nLanes){
    R1_file<-files[grepl("R1",basename(files))][lane]
    R2_file<-files[grepl("R2",basename(files))][lane]
    tempL[[lane]]<-data.frame(
      pool=pool,
      lane=lane,
      rep=1,
      anno=metadata_file,
      fastq_1=R1_file,
      fastq_2=R2_file
    )
  }
  df<-do.call("rbind",tempL)
  out_samplesheet_csv<-paste0(samplesheet_dir,"/",pool,".csv")
  write.csv(df,out_samplesheet_csv,row.names=F,quote=F)
}

