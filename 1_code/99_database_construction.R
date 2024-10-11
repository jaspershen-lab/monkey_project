library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')


dir.create("3_data_analysis/99_database_construction")

setwd("3_data_analysis/99_database_construction")

library(AnnotationForge)

packageVersion("AnnotationForge")

# assignInNamespace(
#   x = ".getNCBIFile",
#   value = function(destDir, fileName, urlStem) {
#     # Use HTTPS instead of FTP
#     url <- paste0("https://ftp.ncbi.nlm.nih.gov/gene/DATA/", fileName)
#     destFile <- file.path(destDir, fileName)
#     if (!file.exists(destFile)) {
#       download.file(url, destFile, mode = "wb")
#     }
#     return(destFile)
#   },
#   ns = "AnnotationForge"
# )

# makeOrgPackageFromNCBI(
#   version = "0.1",
#   author = "xiaotao shen",
#   maintainer = "xiaotao shen",
#   outputDir = ".",
#   tax_id = "9545",
#   genus = "Macaca",
#   species = "nemestrina",
#   rebuildCache = TRUE
# )



# Load the data files
finchFile <- system.file("extdata", "finch_info.txt", package = "AnnotationForge")
finch <- read.table(finchFile, sep = "\t", header = TRUE)

# Prepare gene symbol data
fSym <- finch[, c("GID", "SYMBOL", "GENENAME")]
fSym <- fSym[fSym$SYMBOL != "-", ]
fSym <- fSym[fSym$GENENAME != "-", ]
colnames(fSym) <- c("GID", "SYMBOL", "GENENAME")

# Prepare chromosome data
fChr <- finch[, c("GID", "CHROMOSOME")]
fChr <- fChr[fChr$CHROMOSOME != "-", ]
colnames(fChr) <- c("GID", "CHROMOSOME")

# Load the GO data
finchGOFile <- system.file("extdata", "GO_finch.txt", package = "AnnotationForge")
fGO <- read.table(finchGOFile, sep = "\t", header = TRUE)
fGO <- fGO[fGO$GO != "", ]
fGO <- fGO[fGO$EVIDENCE != "", ]
colnames(fGO) <- c("GID", "GO", "EVIDENCE")
