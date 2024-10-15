library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
data <-
  readr::read_csv(
    "2_data/Results Tables (DIA-NN, library-free search)/PTSamples_Proteograph_Protein_Group_Panel_matrix.csv"
  )

load("3_data_analysis/1_phenotype_data_preparation/phenotype_info.RData")

dir.create("3_data_analysis/2_proteomics_data_preparation",
           showWarnings = FALSE)
setwd("3_data_analysis/2_proteomics_data_preparation/")

variable_info <-
  data %>%
  dplyr::select(Protein.Group:Gene.Names) %>%
  dplyr::mutate(variable_id = paste0("protein_", row_number())) %>%
  dplyr::rename(
    uniprot_id_group = Protein.Group,
    protein_name_group = Protein.Names,
    gene_name_group = Gene.Names
  ) %>%
  dplyr::select(variable_id, everything())

variable_info$uniprot_id <-
  variable_info$uniprot_id_group %>%
  stringr::str_split(";") %>%
  purrr::map_chr(1)

variable_info$protein_name <-
  variable_info$protein_name_group %>%
  stringr::str_split(";") %>%
  purrr::map_chr(1)

variable_info$gene_name <-
  variable_info$gene_name_group %>%
  stringr::str_split(";") %>%
  purrr::map_chr(1)

variable_info$uniprot_id

expression_data <-
  data %>%
  dplyr::select(-c(Protein.Group:Gene.Names))

sample_info <-
  data.frame(sample_id = colnames(expression_data))

length(sample_info$sample_id)

length(phenotype_info$sample_id)

phenotype_info$sample_id <-
  paste0("PT", phenotype_info$sample_id)

intersect(sample_info$sample_id, phenotype_info$sample_id)

setdiff(sample_info$sample_id, phenotype_info$sample_id)
setdiff(phenotype_info$sample_id, sample_info$sample_id)

sample_info$sample_id[sample_info$sample_id == "PTKC7064"] <- "PTK07064"
sample_info$sample_id[sample_info$sample_id == "PTMO6211"] <- "PTM06211"

setdiff(sample_info$sample_id, phenotype_info$sample_id)
setdiff(phenotype_info$sample_id, sample_info$sample_id)

sample_info <-
  sample_info  %>%
  dplyr::left_join(phenotype_info, by = "sample_id")


library(tidymass)

dim(sample_info)
dim(expression_data)
dim(variable_info)

sample_info <-
  sample_info  %>%
  dplyr::mutate(class = case_when(stringr::str_detect(sample_id, "QC") ~ "QC", TRUE ~ "Subject"))

dim(sample_info)
dim(expression_data)

colnames(expression_data) <- sample_info$sample_id

rownames(expression_data) <- variable_info$variable_id

proteomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

###protein number for each sample
dim(proteomics_data)
temp_data <-
  data.frame(
    sample_id = colnames(proteomics_data),
    protein_number = apply(proteomics_data, 2, function(x) {
      sum(!is.na(x))
    })
  ) %>%
  dplyr::arrange(protein_number) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = sample_id))

median(temp_data$protein_number)
mean(temp_data$protein_number)

plot <-
  temp_data %>%
  ggplot(aes(sample_id, protein_number)) +
  geom_segment(aes(xend = sample_id, yend = 0), color = "grey") +
  my_theme +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_hline(yintercept = median(temp_data$protein_number),
             color = "red")  +
  geom_text(aes(
    x = 10,
    y = median(temp_data$protein_number) + 100,
    label = median(temp_data$protein_number)
  ), color = "red")

ggsave(plot,
       filename = "protein_number_per_sample.pdf",
       width = 10,
       height = 7)

###data cleaning
##missing value distribution

plot <-
  massqc::show_missing_values(
    object = proteomics_data,
    show_row_names = FALSE,
    show_column_names = FALSE
  )

plot

###save this plot as pdf
###missing values for each samples
plot <-
  show_sample_missing_values(
    object = proteomics_data,
    percentage = TRUE,
    show_x_text = FALSE,
    show_x_ticks = FALSE,
    desc = TRUE
  ) +
  geom_hline(yintercept = 50, color = "red")
plot
ggsave(plot,
       filename = "missing_values_per_sample.pdf",
       width = 10,
       height = 7)


####missing values for each variable
plot <-
  show_variable_missing_values(
    object = proteomics_data,
    percentage = TRUE,
    show_x_text = FALSE,
    show_x_ticks = FALSE,
    desc = TRUE
  ) +
  geom_hline(yintercept = 50, color = "red")

plot

ggsave(
  plot = plot,
  filename = "missing_values_per_variable.pdf",
  width = 10,
  height = 7
)


#####remove variables that have more than 50% missing values
proteomics_data <-
  proteomics_data %>%
  mutate_variable_na_freq() %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.5)

###if there are outlier samples
show_sample_missing_values(
  object = proteomics_data,
  percentage = TRUE,
  show_x_text = FALSE,
  show_x_ticks = FALSE,
  desc = TRUE
)

####no outlier sampels

###missing value imputation
proteomics_data <-
  impute_mv(object = proteomics_data,
            colnames(proteomics_data),
            method = "knn")

sum(is.na(proteomics_data))

pca_object <-
  proteomics_data %>%
  scale_data() %>%
  run_pca()

plot <-
  pca_score_plot(object = proteomics_data, pca_object = pca_object)

plot

ggsave(plot,
       filename = "pca_score_plot.pdf",
       width = 7,
       height = 7)

save(proteomics_data, file = "proteomics_data.RData")
