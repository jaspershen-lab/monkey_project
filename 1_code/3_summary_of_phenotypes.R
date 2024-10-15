library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
load("3_data_analysis/2_proteomics_data_preparation/proteomics_data.RData")

dir.create("3_data_analysis/3_summary_of_phenotypes", showWarnings = FALSE)
setwd("3_data_analysis/3_summary_of_phenotypes/")

library(tidymass)

sample_info <-
  extract_sample_info(proteomics_data)

sample_info <-
  sample_info  %>%
  dplyr::filter(class == "Subject")

sample_info <-
  sample_info  %>%
  dplyr::filter(!is.na(age))

sample_info$sex
sample_info$age

colnames(sample_info)

plot <-
  sample_info  %>%
  ggplot(aes(x = age)) +
  geom_histogram(aes(fill = sex), color = "black", bins = 100) +
  scale_fill_manual(values = sex_color) +
  my_theme

plot

ggsave(plot,
       filename = "age_sex_distribution.pdf",
       width = 10,
       height = 7)
