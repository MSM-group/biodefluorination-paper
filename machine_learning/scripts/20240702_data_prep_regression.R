temp1 <- read_excel("data/machine_learning/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res3 <- read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "rep3_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res3 <- res3[!is.na(res3)]
