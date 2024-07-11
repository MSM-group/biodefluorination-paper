# Install packages
pacman::p_load("tidyverse", "pROC", "caret", "rsample", "ranger", "microseq", "cowplot",
               "e1071", "dplyr", "seqinr", "readxl", "micropan", "cvms", "DECIPHER")

temp1 <- read_excel("data/machine_learning/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res3 <- read_excel('data/Machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 

expdat <- data.frame(activity = res3) %>%
  bind_cols(label = temp1) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>% # mean activity
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(newlab = paste0("WP_178618037_1_", substr(label, 2, nchar(label))))

# Read in the sequences
m8037 <- read_excel("data/Machine_learning/Batch353c_Robinson_m8037_T7Express.xlsx") %>%
  dplyr::left_join(., expdat, by = c("fd_uname" = "newlab")) %>%
  dplyr::mutate(truth = case_when(activity >= 0.3 ~ "defluor", 
                                  activity <= 0.1 ~ "nondefluor"))  
m9078 <- read_excel("data/Machine_learning/Batch353c_Robinson_m9078_T7Express.xlsx") %>%
  dplyr::mutate(truth = "nondefluor")
comball <- m8037 %>%
  bind_rows(m9078) %>%
  janitor::clean_names() %>%
  dplyr::mutate(nuc = as.character(toupper(dna_sequence))) %>%
  dplyr::mutate(protein = microseq::translate(nuc)) %>%
  dplyr::filter(!is.na(truth))
