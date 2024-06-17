library(readr)
library(stringr)
library(pheatmap)
library(dplyr)
library(tidyverse)

het.inputfile <- read_tsv(paste("data/voila_het.tsv",sep = ""),comment = '#',skip_empty_rows = T)
head(het.inputfile)
dim(het.inputfile)
colnames(het.inputfile)
het <- het.inputfile %>% select(-c(colnames(het.inputfile)[grep("percentile|quantile",colnames(het.inputfile))]))
View(het.inputfile)
colnames(het.inputfile)
#changing_columns <- colnames(het.inputfile)[grep(" changing",colnames(het.inputfile))]

#is.na(het) <- "NA"

#filtro per gli lsv_id con almeno un True tra i 3 test d'ipotesi di changing nei diversi confronti
het <- het.inputfile %>% select(-contains(c("percentile","quantile")),-c("seqid","strand")) %>% 
  filter_at(vars(contains(" changing")), any_vars(str_detect(., pattern ="True")))
View(het)

#het %>% select(matches(" changing"))

