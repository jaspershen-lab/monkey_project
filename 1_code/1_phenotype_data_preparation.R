library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

phenotype_info <-
  readxl::read_xlsx("2_data/Copy of Pigtail blood plasma samples 2023-updated.xlsx")

dir.create("3_data_analysis/1_phenotype_data_preparation",
           showWarnings = FALSE)
setwd("3_data_analysis/1_phenotype_data_preparation/")

phenotype_info$date_collection <-
  as.POSIXlt(lubridate::mdy("03/31/2023"), tz = "UTC")


phenotype_info$age <-
  (phenotype_info$date_collection - phenotype_info$Birth) / 365 %>%
  as.numeric()

phenotype_info <-
  phenotype_info %>%
  dplyr::rename(
    sample_id = "Id",
    cage = Cage,
    species = Species,
    sex = Sex,
    dam = Dam,
    sire = Sire,
    date_birth = Birth,
    group = Group,
    previous_group = "Prev Group",
    box = "Box #",
    labKey_box_whole_blood = "LabKey Box Whole Blood",
    note = "...13",
  )

phenotype_info


