# Install packages
pacman::p_load("tidyverse", "pROC", "caret", "rsample", "ranger", "microseq", "cowplot", "ips",
               "e1071", "dplyr", "seqinr", "readxl", "micropan", "cvms", "DECIPHER")
library("ips")

# Read in the hads to classify
had3 <- readAAStringSet("data/machine_learning/50_unique_validation_hads.fasta")
had3al <- AlignSeqs(had3)
width(had3al)
BrowseSeqs(had3al)
# Try aligning profiles
alprof <- readAAStringSet("data/Machine_learning/mutlib_all_seqs_aligned.fasta")
BrowseSeqs(alprof)
width(alprof)
aaal <- AlignSeqs(AAStringSet(c(alprof, had3al)))
aaal2 <- AlignProfiles(alprof, had3al[2])
BrowseSeqs(aaal2)

# Read in the new hads to classify 
guop <- readAAStringSet("data/machine_learning/validation/guopingia_1000_hits.txt")
duplicated(guop)
gordoni <- readAAStringSet("data/machine_learning/validation/gordoni_1000_hits.txt")
accs <- word(names(guop), sep = " ", 1)
writeLines(accs, "data/machine_learning/validation/guopingia_1000_accs.txt")
accs2 <- word(names(gordoni), sep = " ", 1)
inter <- intersect(accs, accs2)
writeLines(accs2, "data/machine_learning/validation/gordoni_1000_accs.txt")

# Read in rawdata
dat2 <- read_csv("data/machine_learning/20240701_defluorinases_entire_alignment_for_classification.csv") %>%
  dplyr::slice(1:27)# sample testing set
dat3 <- dat2 %>%
  dplyr::slice(1:27) %>%
  dplyr::select(-nams, -truth)
dim(dat3)

modelr <- readRDS("data/Machine_learning/20240701_classification_random_forest.rds")
modelr
dat2$truth

# Make a prediction
modelr$
newpred <- predict(modelr, newdata = dat3)
length(newpred)
truthvec <- as.factor(dat2$truth)
truthvec
newpred
confusionMatrix(newpred, truthvec)

