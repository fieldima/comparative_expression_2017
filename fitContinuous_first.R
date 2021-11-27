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
adequacy <- function (tree, part) {
  result <- vector("list", length = length(tree))
  count <- 1
  for (tr in tree){
    phylo <- tr@phylo
    dat <- tr@data %>% select(label, Averaged.RPKM.brain:Averaged.RPKM.liver) %>%
      slice(1:((nrow(tr@data)+1)/2)) %>% column_to_rownames(var = "label") %>%
      as.matrix()
    dat[is.na(dat)] = 0
    fitBM <- fitContinuous(phylo, dat = dat[,part], model = "BM")$opt
    fitOU <- fitContinuous(phylo, dat = dat[,part], model = "OU")$opt
    fitEB <- fitContinuous(phylo, dat = dat[,part], model = "EB")$opt
    aic <- c(fitBM[["aic"]], fitOU[["aic"]], fitEB[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, "BM")), 
            ifelse(min(aic) == aic[2], list(c(fitOU, "OU")), 
            list(c(fitEB, "EB"))))
    result[[count]] <- fit
    count <- count + 1
  }
  result
}

first_2k <- calibrated[1:2000]

#Measure adequacy
brain <- adequacy(first_2k, 1)
saveRDS(brain, file = "arbutus/brain_fit")

heart <- adequacy(first_2k, 2)
saveRDS(heart, file = "arbutus/heart_fit")

kidney <- adequacy(first_2k, 3)
saveRDS(kidney, file = "arbutus/kidney_fit")

testis <- adequacy(first_2k, 4)
saveRDS(testis, file = "arbutus/testis_fit")

cerebellum <- adequacy(first_2k, 5)
saveRDS(cerebellum, file = "arbutus/cerebellum_fit")

liver <- adequacy(first_2k, 6)
saveRDS(liver, file = "arbutus/liver_fit")


save.image(file = "arbutus/fitContinuous.RData")
