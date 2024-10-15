library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
load("3_data_analysis/2_proteomics_data_preparation/proteomics_data.RData")

dir.create("3_data_analysis/4_markers_of_ages", showWarnings = FALSE)
setwd("3_data_analysis/4_markers_of_ages/")

library(tidymass)

expression_data <-
  extract_expression_data(proteomics_data)

sample_info <-
  extract_sample_info(proteomics_data)

variable_info <-
  extract_variable_info(proteomics_data)

####use the linear mixed model to find the proteins that are associated with age,
###and sex should be included as a fixed effect

plot_age_male <-
  sample_info %>%
  dplyr::filter(sex == "Male") %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  ggplot(aes(sample_id, age)) +
  geom_point() +
  my_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(plot_age_male,
       filename = "plot_age_male.pdf",
       width = 10,
       height = 7)

plot_age_female <-
  sample_info %>%
  dplyr::filter(sex == "Female") %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = sample_id)) %>%
  ggplot(aes(sample_id, age)) +
  geom_point() +
  my_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot_age_female
ggsave(
  plot_age_female,
  filename = "plot_age_female.pdf",
  width = 10,
  height = 7
)

# age_markers <-
#   vector("list", length = nrow(variable_info))
#
# for (i in 1:nrow(variable_info)) {
#   cat(i, " ")
#   temp_data <-
#     data.frame(
#       subject_id = sample_info$sample_id,
#       age = sample_info$age,
#       sex = sample_info$sex,
#       value = as.numeric(expression_data[i, ])
#     )
#
#   temp_data_male <-
#     temp_data %>%
#     dplyr::filter(sex == "Male")
#
#   temp_data_female <-
#     temp_data %>%
#     dplyr::filter(sex == "Female")
#
#   cor_test_male <- cor.test(temp_data_male$age, temp_data_male$value)
#   cor_test_female <- cor.test(temp_data_female$age, temp_data_female$value)
#
#   age_markers[[i]] <-
#     data.frame(
#       variable_id = variable_info$variable[i],
#       cor_male = cor_test_male$estimate,
#       p_male = cor_test_male$p.value,
#       cor_female = cor_test_female$estimate,
#       p_female = cor_test_female$p.value
#     )
#
# }
#
#
# save(age_markers, file = "age_markers.RData")
load("age_markers.RData")

age_markers <-
  age_markers %>%
  dplyr::bind_rows()

rownames(age_markers) <- NULL

age_markers$p_female_fdr <- p.adjust(age_markers$p_female, method = "fdr")
age_markers$p_male_fdr <- p.adjust(age_markers$p_male, method = "fdr")

library(ggside)

test <-
  cor.test(age_markers$cor_male, age_markers$cor_female)

plot <-
  age_markers %>%
  ggplot(aes(cor_male, cor_female)) +
  geom_abline(color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(alpha = 0.5) +
  ggside::geom_xsidepoint(aes(color = cor_male, size = -log(p_male_fdr, 10))) +
  ggside::geom_ysidepoint(aes(color = cor_female, size = -log(p_female_fdr, 10))) +
  scale_color_gradient2(low = "blue",
                        mid = "white",
                        high = "red") +
  scale_size_continuous(range = c(0.1, 3)) +
  my_theme +
  labs(x = "Correlation (Male)", y = "Correlation (Female)") +
  geom_text(
    x = -0.5,
    y = 0.4,
    label = paste(
      "correlation",
      round(test$estimate, 4),
      "\np-value:",
      signif(test$p.value, 2)
    ),
    color = "red"
  )

plot

ggsave(plot,
       filename = "age_markers.pdf",
       width = 9.5,
       height = 7)


####markers that are consistent in male and female
age_markers_consistent_up <-
  age_markers %>%
  dplyr::filter(cor_female > 0 & p_female_fdr < 0.05 &
                  cor_male > 0 & p_male_fdr < 0.05)


age_markers_consistent_down <-
  age_markers %>%
  dplyr::filter(cor_female < 0 & p_female_fdr < 0.05 &
                  cor_male < 0 & p_male_fdr < 0.05)

dim(age_markers_consistent_up)
dim(age_markers_consistent_down)



