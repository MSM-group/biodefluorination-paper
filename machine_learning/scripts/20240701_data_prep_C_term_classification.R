# Install packages
pacman::p_load("tidyverse", "pROC", "caret", "rsample", "ranger", 
               "ggseqlogo", "microseqlogo", "seqinr", "readxl", "DECIPHER")

# Read in the template
temp1 <- read_excel("data/Machine_learning/template_linearized.xlsx", col_names = F) %>%
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

res3 <- res3[!is.na(res3)]

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

mutlib <- AAStringSet(x = comball$protein)
names(mutlib) <- comball$fd_uname

# Align all sequences
pos <- readAAStringSet("data/Machine_learning/pos_defluorinases_deduplicated.fasta")
neg <- readAAStringSet("data/Machine_learning/neg_defluorinases_deduplicated.fasta")
comb <- AlignSeqs(AAStringSet(c(pos, neg, mutlib)))
names(comb) <- c(names(pos), names(neg), names(mutlib))
# BrowseSeqs(comb)
# writeXStringSet(comb, "data/Machine_learning/C_term_all_seqs_aligned.fasta")

# Load alignment file
seqs <- read.alignment("data/Machine_learning/mutlib_all_seqs_aligned.fasta", format = "fasta") 
lb_rham <- seqs$seq[[1]] 

# Identify the start and end of the C-terminal region
nchar("PTLYVPRPLEYGAVNNHVEPEAKYDHETVADFRELAARLGV") # full region is 41 amino acids, however alignment has a gap
tri.pos1 <- words.pos("ptl", lb_rham) # C-terminal motif starts with 'PTL'
end.pos1 <- words.pos("nnh", lb_rham)
tri.pos1
end.pos1

# Section 1
nuc1 <- lapply(seqs$seq, function(x) { substr(x,tri.pos1,end.pos1) })
nuc1
nucr1 <- unlist(nuc1)
nucr1

# Section 2
tri.pos2 <- words.pos("vepe",lb_rham) # ptl #mdvsnv
end.pos2 <- words.pos("gv-", lb_rham)
tri.pos2
end.pos2
nuc2<-lapply(seqs$seq,function(x) { substr(x,tri.pos2,end.pos2+1) })
nuc2
nucr2<-unlist(nuc2)

# Concatenate two sections
nucr <- paste0(nucr1, nucr2)
names(nucr) <- seqs$nam
appendf <- data.frame(nams = names(nucr), motif = nucr)

# 206 - 221
# EYGAVNNHHEVPEAKY
# Make a logo of the C-terminal motif 
motif_seqs <- data.frame(toupper(substr(appendf$motif, 
                                start = 23,
                                stop = nchar(appendf$motif))))

motif_seqs$toupper.substr.appendf.motif..start...23..stop...nchar.appendf.motif...
fig1 <- ggplot() + 
  geom_logo(motif_seqs, method = "p", 
            col_scheme = "auto") +
  theme_logo() +
  theme(axis.text.x = element_text(angle = 45))
fig1
fig1$scales$scales[[1]] <- scale_x_continuous(breaks = c(1:18), labels=as.character(c(274:(275+16))))
fig1

ggsave(fig1, file = "output/shortened_C_terminal_motif.png", width = 6, height = 3)

# Merge with the defluorination data
merg <- appendf %>%
  dplyr::mutate(motif = toupper(motif))
nchar(merg$motif)
tri.pos1
end.pos1

tri.pos2
end.pos2
merg_split <- merg %>%
  tidyr::separate(motif, into = paste0("residue", c((tri.pos1-12):(end.pos1-11), (tri.pos2-11):(290-11), "_", (292-11):(end.pos2-10))), sep = "")
  #tidyr::separate(motif, into = paste0("residue", 245:(245 + width(nucr)[1])), sep = "") 
  #tidyr::separate(motif, into = paste0("residue", 197:(197 + width(nucr)[1])), sep = "") 

appendf <- data.frame(nams = names(nucr), motif = nucr)


# Write unencoded data frame to file
reg_df <- merg_split %>%
  dplyr::mutate(truth = c(rep("defluor", length(pos)), rep("nondefluor", length(neg)),
                          comball$truth))
write_csv(reg_df, "data/machine_learning/20240701_defluorinases_C_term_for_classification.csv")

# Remove columns that the machine learning model should not learn e.g., nams, truth
rawdat <- reg_df %>%
  select(-(1:2)) 
colnames(rawdat) <- paste0(colnames(rawdat), "_", as.character(rawdat[1,]))
colnames(rawdat)[grepl("truth", colnames(rawdat))] <- "truth"
colnames(rawdat)[grepl("nams", colnames(rawdat))] <- "nams"

# Remove variables with nonzero variance 
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 

dat <- rawdat  %>%
  select(-all_of(which_rem))

# Check and remove duplicates
dat <- dat[!duplicated(dat),] %>%
  dplyr::filter(!is.na(truth))
rownames(dat)
colnames(dat)
# Set random seed 
set.seed(20240521) # This is the one

# Split into test and training data - random option
dat_split <- rsample::initial_split(dat, prop = 0.8, 
                                    strata = "truth")
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "truth")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "truth")]

# Dependent variable
y_train <- dat_train$truth
y_test <- dat_test$truth

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$nams)

# Optional tuning of random forest parameters
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0), ncol(df_train))
mtrys # number of variables available for splitting at each tree node

# Tuning parameter grid
rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("gini", "extratrees"),
                       min.node.size = 1)

# Train with 10-fold cross validation and 3 repetitions
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 3,  
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
rf$finalModel$prediction.error # Out-of-bag error
saveRDS(rf, "data/Machine_learning/20240701_C_term_classification_random_forest.rds")

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp

pdf("data/Machine_learning/C_terminal_only_importance_plot_full_length_classification_top20.pdf", width = 5, height = 3)
rf1 <- ggplot(rf_imp, top = 20) + 
  xlab("") +
  theme_classic()
rf1
dev.off()

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred
see_preds <- bind_cols(rownames(dat_test), form_test$y_test, rf_pred)
see_preds
cm_rf <- confusionMatrix(rf_pred, as.factor(y_test))

# ROC curve
rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                    predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                    ci = T)

pl1 <- plot(rf_roc, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.95 
     , print.auc.y = 0.8,
     xlim = c(1.1,-0.1), ylim = c( 0, 1.1))
pl1

rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                   predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                   ci = F)

pdf("data/Machine_learning/C_term_only_auroc_curve_classification.pdf", width = 3, height = 3)
plot(rf_roc, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.95 
     , print.auc.y = 0.8,
     xlim = c(1.1,-0.1), ylim = c( 0, 1.1))
dev.off()



