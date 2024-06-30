# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", 
               "colorspace", "cowplot", "tidymodels", "ranger", "tree", 
               "rsample", "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Read in the dataset
reg_df <- read_csv('data/machine_learning/20240629_defluorinases_entire_alignment_for_classification.csv')

rawdat <- reg_df %>%
  select(-2)

# Remove variables with nonzero variance (optional)
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, uniqueCut = 0.999999)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
length(which_rem) # cut 73

# Check for duplicates
dat <- rawdat  %>%
  select(-which_rem)
dat <- dat[!duplicated(dat),] # 107 rows
nrow(dat)
colnames(dat) <- paste0(colnames(dat), "_", dat[1,])
colnames(dat)[ncol(dat)] <- "truth"
colnames(dat)[1] <- "nams"

# Read in the random forest model from cross-validation
rf_20200629 <- readRDS("data/machine_learning/20240629_classification_random_forest.rds")
rf_20200629$bestTune$mtry
rf_20200629$bestTune$splitrule
rf_20200629$bestTune$min.node.size

# Set random seed 
set.seed(5678)  

rf_models <- vector(mode = "list",
                    length = 1000) # 1000 different splits

# Split randomly into 10 test/training sets
for(i in 1:1000) { 
  # Split into test and training data
  dat_split <- initial_split(dat, strata = "truth")
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-nams, -truth)
  x_test <- dat_test %>%
    dplyr::select(-nams, -truth)
  
  y_train <- as.factor(dat_train$truth)
  y_test <- as.factor(dat_test$truth)
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$nams)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)
  
  ## Train model using set parameters (determined by cross-validation)
  rf_1 <- ranger(y_train ~., data = form_train, 
                 num.trees = 1000, 
                 splitrule = "extratrees", # rf_20200629$bestTune$splitrule
                 mtry = 8, # rf_20200629$bestTune$mtry
                 min.node.size = 1, # rf_20200629$bestTune$min.node.size,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classifcation accuracy 1 - OOB error

  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$confmat <- confusionMatrix(rf_1_pred$predictions, as.factor(dat_test$truth))
}

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$oob }) 
mean_training_vec
traindf <- tibble(mean_training_vec)
mean(mean_training_vec)
sd(mean_training_vec) # 5.74%
write_tsv(traindf,"data/machine_learning/20240629_classification_mean_training_accuracy_1000_test_train_splits.txt")

p1 <- ggplot(traindf) +
  geom_histogram(aes( x = mean_training_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean training set accuracy") +
  ylab("Frequency")
p1

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
testdf <- tibble(mean_testing_vec)

p2 <- ggplot(testdf) +
  geom_histogram(aes( x = mean_testing_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean testing set accuracy") +
  ylab("Frequency")
p2

write_tsv(tibble(mean_testing_vec),"data/machine_learning/20240629_classification_mean_testing_accuracy_1000_test_train_splits.txt")
mean(mean_training_vec)
sd(mean_training_vec) 

# Combine into a single histogram
traindf$split <- "training"
colnames(traindf) <- c("score", "split")
trainsumm <- data.frame(tibble(t(summary(traindf$score))), stringsAsFactors = F)
trainsumm

testdf$split <- "testing"
colnames(testdf) <- c("score", "split")
summary(testdf$score)
testsumm <- data.frame(tibble(t(summary(testdf$score))), stringsAsFactors = F)
testsumm


allsumm <- data.frame(c(trainsumm,testsumm)) 
write_csv(allsumm, "data/machine_learning/20240629_classification_train_test_1000_summaries.csv")


alldf <- traindf %>%
  bind_rows(testdf)
head(alldf)
alldf$split <- as.factor(alldf$split)
alldf$score
pdf("data/machine_learning/20240629_1000_binary_classification_training_testing_split_hist.pdf", width = 7, height = 4)
p3 <- ggplot(alldf, aes(x= score, fill = split)) +
  geom_histogram(alpha = 0.5, binwidth = 0.01,
                  position = 'identity') +
  theme_pubr(base_size = 12) +
  xlab("Classification accuracy") +
  ylab("Count") #+
 # scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
p3
dev.off()

mean_F1_vec <- sapply(1:length(rf_models), function(x) {  rf_models[[x]]$confmat$byClass[['F1']] }) 
mean_F1_vec
dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) %>%
  bind_cols(tibble(mean_F1_vec))
dat_df
write_csv(dat_df, "data/Machine_learning/20240629_full_alignment_classification_accuracy_1000_training_test_combined.csv")


sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
colMeans(dat_df)
sd_df

rf_models[[1]]$confmat

