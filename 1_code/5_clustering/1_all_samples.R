library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
load("3_data_analysis/2_proteomics_data_preparation/proteomics_data.RData")

dir.create(
  "3_data_analysis/5_clustring/1_all_samples",
  showWarnings = FALSE,
  recursive = TRUE
)
setwd("3_data_analysis/5_clustring/1_all_samples/")

library(tidymass)

proteomics_data <-
  proteomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(!is.na(age))

expression_data <-
  extract_expression_data(proteomics_data)

sample_info <-
  extract_sample_info(proteomics_data)

variable_info <-
  extract_variable_info(proteomics_data)


###combine samples into different age groups
sample_info$age_group <-
  cut(sample_info$age,
      breaks = c(0:24),
      labels = c(paste(0:23, "-", 1:24, sep = "")))

expression_data_new <-
  as.character(sort(unique(sample_info$age_group))) %>%
  purrr::map(function(x) {
    cat(x, " ")
    idx <-
      which(sample_info$age_group == x)
    expression_data[, idx, drop = FALSE] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_new) <-
  as.character(sort(unique(sample_info$age_group)))

sample_info_new <-
  data.frame(
    sample_id = colnames(expression_data_new),
    age = stringr::str_split(colnames(expression_data_new), "-") %>%
      lapply(function(x) {
        mean(as.numeric(x))
      }) %>% unlist
  )


###fuzzy c meaning clustering
library(Mfuzz)
time <- sample_info_new$age

temp_data <- rbind(time, expression_data_new)
rownames(temp_data)[1] <- "time"
rownames(temp_data)

write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")
data.s <- standardise(data)
# data.s <- data
m1 <- mestimate(data.s)
m1
#
# plot <-
#   Dmin(
#     data.s,
#     m = m1,
#     crange = seq(2, 40, 2),
#     repeats = 3,
#     visu = TRUE
#   )
#
# plot
#
# plot <-
#   plot %>%
#   data.frame(distance = plot, k = seq(2, 40, 2)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21,
#              size = 4,
#              fill = "black") +
#   # geom_smooth() +
#   geom_segment(aes(
#     x = k,
#     y = 0,
#     xend = k,
#     yend = distance
#   )) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(x = "Cluster number", y = "Min. centroid distance") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
#
# plot
#
# ggsave(plot,
#        filename = "distance_k_number.pdf",
#        width = 7,
#        height = 7)

cluster_number <- 12

c <- mfuzz(data.s, c = cluster_number, m = m1)

save(c, file = "c")
load("c")

# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
# center <- c$centers

membership_cutoff <- 0.8

center <-
  get_mfuzz_center(data = data.s,
                   c = c,
                   membership_cutoff = 0.5)

rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

mfuzz.plot(
  eset = data.s,
  min.mem = 0.8,
  cl = c,
  mfrow = c(3, 4),
  time.labels = time,
  new.window = FALSE
)

library(ComplexHeatmap)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

variable_info <-
  proteomics_data@variable_info

temp_data <-
  expression_data_new %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

for (idx in 1:cluster_number) {
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    # dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    center[idx, , drop = TRUE] %>%
    unlist() %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = as.numeric(time))
  
  temp <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "sample_id",
      values_to = "value"
    ) %>%
    dplyr::left_join(sample_info_new[, c("sample_id", "age")], by = "sample_id") %>%
    dplyr::mutate(age = as.numeric(age))
  
  plot <-
    temp %>%
    dplyr::arrange(desc(membership)) %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(age, value, group = variable_id)) +
    geom_line(aes(color = membership), alpha = 0.7) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.justification = c(0, 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z - score",
      title = paste("Cluster ", idx, " (", nrow(cluster_data), " proteincs)", sep = "")
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0)
  # viridis::scale_color_viridis()
  
  plot
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
}

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      # dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff)
    temp <-
      temp %>%
      dplyr::mutate(cluster_raw = cluster) %>%
      dplyr::mutate(cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

cluster_info %>%
  dplyr::count(cluster)

cluster_info %>%
  dplyr::filter(membership > 0.8) %>%
  dplyr::count(cluster)

final_cluster_info <-
  cluster_info

save(final_cluster_info, file = "final_cluster_info")
write.csv(final_cluster_info, "final_cluster_info.csv", row.names = FALSE)
