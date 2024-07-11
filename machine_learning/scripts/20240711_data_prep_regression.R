# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggpmisc", "caret", 
               "colorspace", "cowplot", "tidymodels", "ranger", "tree", 
               "rsample", "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Read in the template
temp1 <- read_excel("data/Machine_learning/template_linearized.xlsx",
                    col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res <- read_excel('data/Machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                  col_names = F, sheet = "average_linearized") %>%
                   #col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 

res <- res[!is.na(res)]

# Delta calculation relative to WT
wt <- res[temp1 == "WT"]
wt

# Delta
delta <- wt - res # delta between WT -> Ala_mut
dat <- bind_cols(label = temp1, 
                 delta = delta) %>%
  dplyr::filter(label != "WT") %>%
  dplyr::mutate(position_index = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::mutate(position_from = substr(label, 2, 2)) %>%
  dplyr::mutate(position_to = substr(label, nchar(label), nchar(label))) %>%
  dplyr::mutate(unique_id = paste0(position_from, "_", position_index, "_", position_to, "_fwd"))

props <- read_csv("data/5_aa_properties.csv")
colnames(props) <- word(colnames(props), sep = "_", -1)
colnames(props)[1] <- "AA_ABREV"
colnames(props)
feat_df <- dat %>%
  dplyr::left_join(., props, by = c("position_from" = "AA_ABREV")) %>%
  dplyr::mutate(label = rownames(.)) %>%
  dplyr::select(-position_from, -position_to)
feat_df

# Try predicting the activity of 20%
set.seed(78910)

# Split into test and training data 
# 80% training
# 20% test
dat <- feat_df 
dat_split <- rsample::initial_split(dat, prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) 

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("unique_id", "label", "delta")]
x_test <- dat_test[,!colnames(dat_test) %in% c("unique_id", "label", "delta")]
colnames(x_train)

# Dependent variable
y_train <- dat_train$delta
y_test <- dat_test$delta
y_test # check there is a mix of your two variables

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$unique_id)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$label)

# Optional tuning of random forest parameters
mtrys <- c(round(log2(ncol(df_train)), 0), round(ncol(df_train), 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("extratrees"),
                       min.node.size = 1)

# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 3, # increase this to 3 when you run the code 
                           verboseIter = T, 
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
rf$finalModel$prediction.error # out-of-bag error


#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp

#pdf("data/machine_learning/variable_importance_plot.pdf", width = 10)
ggplot(rf_imp, top = 5) + 
  xlab("") +
  theme_pubr(base_size = 20)

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$delta))
cm_rf

rfdf_train <- data.frame(cbind(rf$finalModel$predictions, y_train), stringsAsFactors = F)
colnames(rfdf_train) <- c("Pred", "Obs")


rfdf_test <- data.frame(cbind(rf$finalModel$predictions, y_train), stringsAsFactors = F)
colnames(rfdf_test) <- c("Pred", "Obs")

ggplot(rfdf_train, aes(x = Pred, y = Obs)) +
  geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  # geom_smooth(se = FALSE, col = "red", method = "lm",
  #              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = TRUE)

ggplot(rfdf_test, aes(x = Pred, y = Obs)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  # geom_smooth(se = FALSE, col = "red", method = "lm",
  #             lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  xlim(-0.4, 0.4) +
  ylim(-0.4, 0.4)



