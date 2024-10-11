library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
data <-
  readr::read_csv(
    "2_data/Results Tables (DIA-NN, library-free search)/PTSamples_Proteograph_Protein_Group_Panel_matrix.csv"
  )

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


library(biomaRt)


# Load the packages
library(AnnotationHub)
library(AnnotationForge)

# Create an AnnotationHub object
hub <- AnnotationHub()

# Search for Macaca nemestrina data
q <- AnnotationHub::query(hub, "Macaca nemestrina")

print(q)

id <- q$ah_id[length(q)]



library(AnnotationForge)

makeOrgPackageFromNCBI(
  version = "0.1",
  author = "xiaotao shen",
  maintainer = "xiaotao shen",
  outputDir = ".",
  tax_id = "9545",
  genus = "Macaca",
  species = "nemestrina",
  rebuildCache = TRUE
)


###our samples are from Macaca nemestrina (pig-tailed macaque)

library(AnnotationHub)
hub <- AnnotationHub()
q <- AnnotationHub::query(x = hub, pattern = "Macaca nemestrina")
id <- q$ah_id[length(q)]
macaca_nemestrina <- hub[[id]]

library(clusterProfiler)

bitr(
  geneID = variable_info$uniprot_id[1],
  fromType = "UNIPROTID1",
  toType = "ENTREZID",
  OrgDb = macaca_nemestrina
)

expression_data <-
  data %>%
  dplyr::select(-c(Protein.Group:Gene.Names))
