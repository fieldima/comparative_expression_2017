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

toy <- calibrated[1:4]

#Tests
test1 <- adequacy(toy, 1)
test2 <- adequacy(toy, 2)
test3 <- adequacy(toy, 3)
test4 <- adequacy(toy, 4)
test5 <- adequacy(toy, 5)
test6 <- adequacy(toy, 6)


#Measure adequacy
brain <- adequacy(calibrated, 1)
heart <- adequacy(calibrated, 2)
kidney <- adequacy(calibrated, 3)
testis <- adequacy(calibrated, 4)
cerebellum <- adequacy(calibrated, 5)
liver <- adequacy(calibrated, 6)

save.image(file = "arbutus/fitContinuous.RData")
