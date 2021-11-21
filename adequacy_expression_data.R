#Arbutus on comparative expression data

#Load packages
library(geiger)
library(arbutus)
library(treeio)
library(ggtree)
library(tidyverse)

#Load treedata
calibrated <- readRDS("export/calibrated_trees")

#Make a function to iterate through each list
adequacy <- function (tree) {
  for (tr in tree){
    tib <- as_tibble(tr)
    phylo <- tib %>% select(parent, node, branch.length, label.x) %>%
      mutate(label = label.x) %>% select(-label.x)
    dat <- tib %>% select(label.x, Averaged.RPKM.brain:Averaged.RPKM.liver)
  }
}

#Pull entry from list
cal_entry <- calibrated[[1]]

#Convert to tibbles
cal_tib <- as_tibble(cal_entry)

#Extract phylo objects
cal_phylo <- cal_tib %>% select(parent, node, branch.length, label.x) %>%
  mutate(label = label.x) %>% select(-label.x) %>% as.phylo()

#Extract data
cal_data <- cal_tib %>% select(label.x, Averaged.RPKM.brain:Averaged.RPKM.liver)
cal_dat <- cal_data %>% slice(1:((nrow(cal_data)+1)/2)) %>% column_to_rownames(var = "label.x") %>%
  as.matrix()

#Run arbutus
cal_results <- arbutus(cal_phylo, data = cal_dat[,"Averaged.RPKM.heart"])
