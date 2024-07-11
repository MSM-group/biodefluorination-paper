# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggpmisc", "caret", "gplots",
               "colorspace", "cowplot", "tidymodels", "ranger", "tree",
               "rsample", "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Read in the template
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

res4 <- read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "rep4_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res4 <- as.numeric(res4[!is.na(res4)])

# Delta calculation relative to WT
wt3 <- res3[temp1 == "WT"]
wt3

wt4 <- as.numeric(res4[temp1 == "WT"])
res4

# Delta
delta3 <- wt3 - res3 # delta between WT -> Ala_mut
delta4 <- wt4 - res4

delta <- c(delta3, delta4)
label <- c(temp1, temp1) # check the NAs

# Delta
dat <- bind_cols(label = label, 
                 delta = delta) %>%
  dplyr::filter(label != "WT") %>%
  dplyr::mutate(position_index = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::mutate(position_from = substr(label, 2, 2)) %>%
  dplyr::mutate(position_to = substr(label, nchar(label), nchar(label))) %>%
  dplyr::mutate(unique_id = paste0(position_from, "_", position_index, "_", position_to, "_fwd")) %>%
  dplyr::mutate(equal_segment = ntile(position_index, n = 5)) %>%
  dplyr::mutate(chimera_segment_six = ntile(position_index, n = 6)) %>%
  dplyr::mutate(chimera_segment = case_when(position_index %in% 1:35 ~ 1,
                                            position_index %in% 36:69 ~ 2,
                                            position_index %in% 70:109 ~ 3,
                                            position_index %in% 110:150 ~ 4,
                                            position_index %in% 151:196 ~ 5,
                                            position_index %in% 197:238 ~ 6,
                                            TRUE ~ NA))

p1 <-  ggplot(dat, aes(x = position_index, y = delta)) +
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment),alpha = 0.95)) +
  theme_pubr() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  xlab("Protein position") +
  ylab("Enzyme activity: wild-type - variant")+
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) 
p1  

p2 <- ggplot(dat, aes(x = as.factor(chimera_segment), y = delta)) +
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) + 
  geom_violin(aes(group = as.factor(chimera_segment), color = as.factor(chimera_segment), fill = as.factor(chimera_segment)), alpha = 0.2) +
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Chimera segment") +
  ylab("Enzyme activity: wild-type - variant")+
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) +
  geom_signif(comparisons = list(c(6,5),
                                 c(6,4),
                                 c(6,3),
                                 c(6,2),
                                 c(6,1)),
              y_position = 0.8, tip_length = 0, vjust = 0.7,
              step_increase = 0.05,  # To stagger the bars
              textsize = 5,
              map_signif_level = T)
p2



ggsave(p1, file = "output/distribution_plot.png", height = 4, width = 3)
ggsave(p2, file = "output/violin_plot.png", height = 6, width = 3)

