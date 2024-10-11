library(tidyverse)
library(ggplot2)

library(httr)
library(jsonlite)

map_uniprot_ids <- function(uniprot_ids, to_db) {
  ids <- paste(uniprot_ids, collapse = ",")
  
  # Submit the mapping job to UniProt
  response <- POST(
    url = "https://rest.uniprot.org/idmapping/run",
    body = list(from = "UniProtKB_AC-ID", to = to_db, ids = ids),
    encode = "form"
  )
  
  # Extract the job ID from the response
  job <- content(response)
  job_id <- job$jobId
  
  # Poll the status of the job until it's finished
  repeat {
    Sys.sleep(3)
    status_response <- GET(paste0("https://rest.uniprot.org/idmapping/status/", job_id))
    status_content <- content(status_response)
    status <- status_content$jobStatus
    if (status == "FINISHED")
      break
  }
  
  # Retrieve the results
  result_response <- GET(paste0(
    "https://rest.uniprot.org/idmapping/results/",
    job_id,
    "?format=json"
  ))
  result_content <- content(result_response, as = "text")
  result_json <- fromJSON(result_content)
  
  # Extract mappings into a data frame
  mappings <- result_json$results
  if (length(mappings) > 0) {
    df <- data.frame(
      UniProtID = sapply(mappings, function(x)
        x$from),
      MappedID = sapply(mappings, function(x)
        x$to),
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      UniProtID = uniprot_ids,
      MappedID = NA,
      stringsAsFactors = FALSE
    )
  }
  
  return(df)
}

# # Your UniProt IDs
# uniprot_ids <- c("A0A0K2ZAK1", "A0A0U5NNJ0", "A0A0U5P8Q9")
# 
# # Map to Ensembl Gene IDs
# ensembl_mapping <- map_uniprot_ids(uniprot_ids, to_db = "ENSEMBL_ID")

# 
# 
# library(mygene)
# 
# 
# library(mygene)
# 
# uniprot_ids <- c("A0A0K2ZAK1", "A0A0U5NNJ0", "A0A0U5P8Q9")
# 
# # Query MyGene.info for mappings
# results <- queryMany(
#   uniprot_ids,
#   scopes = "uniprot",
#   fields = c("entrezgene", "ensembl.gene"),
#   species = 9545  # Use Taxonomy ID instead of species name
# )
