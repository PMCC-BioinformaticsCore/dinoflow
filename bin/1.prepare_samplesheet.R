library(R.utils)
library(tidyverse)

#fetch the id to create samplesheets

project_name<-"nsd_pilot"
platform<-"NovaSeq"
run_id<-"230712_A01524_0152_AH5HWKDRX3"
researcher<-"Lizzy Pijpers"

samplesheet_dir<-paste0("/researchers/nenad.bartonicek/projects/mac-seq/sample_sheets/",project_name)
metadata_dir<-paste0("/researchers/nenad.bartonicek/projects/mac-seq/metadata/",project_name)
metadata_file<-paste0(metadata_dir,"/",project_name,".csv")

dfL<-list()

in_dir<-paste0(
  "/pipeline/Archives/",
  platform,
  "/",run_id,
  "/ProjectFolders",
  "/Project_",gsub(" ","-",researcher),
  "/"
)

pool_dirs<-list.files(
  in_dir,full.names=T
)

for(pool_dir in pool_dirs){
  pool_name<-gsub("Sample_","",basename(pool_dir))
  pool_files<-list.files(
    pool_dir,full.names=T
  )
  nLanes<-sum(grepl("R1",basename(pool_files)))
  tempL<-list()
  for(lane in 1:nLanes){
    R1_file<-pool_files[grepl("R1",basename(pool_files))][lane]
    R2_file<-pool_files[grepl("R2",basename(pool_files))][lane]
    tempL[[lane]]<-data.frame(
      pool=pool_name,
      lane=lane,
      rep=1,
      anno=metadata_file,
      fastq_1=R1_file,
      fastq_2=R2_file
    )
  }
  df<-do.call("rbind",tempL)
  out_samplesheet_csv<-paste0(samplesheet_dir,"/",pool_name,".csv")
  write.csv(df,out_samplesheet_csv,row.names=F,quote=F)
}

