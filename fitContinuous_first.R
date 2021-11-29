#Arbutus on comparative expression data

#Load packages
library(geiger)
library(arbutus)
library(treeio)
library(ggtree)
library(tidyverse)

#Load treedata
calibrated <- readRDS("export/calibrated_trees") %>% discard(~ Nnode2(.x) < 9)

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
    fitBM <- fitContinuous(phylo, dat = dat[,part], model = "BM")
    fitOU <- fitContinuous(phylo, dat = dat[,part], model = "OU")
    fitEB <- fitContinuous(phylo, dat = dat[,part], model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
            ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                   list(c(fitEB, model = "EB"))))
    result[[count]] <- fit
    count <- count + 1
  }
  result
}

first <- calibrated[1:2777]
second <- calibrated[2778:5554]
third <- calibrated[5555:8333]

#Measure adequacy
brain <- adequacy(first, 1)
saveRDS(brain, file = "arbutus/brain_fit_first")
brain <- adequacy(second, 1)
saveRDS(brain, file = "arbutus/brain_fit_second")
brain <- adequacy(third, 1)
saveRDS(brain, file = "arbutus/brain_fit_third")


heart <- adequacy(first, 2)
saveRDS(heart, file = "arbutus/heart_fit_first")
heart <- adequacy(second, 2)
saveRDS(heart, file = "arbutus/heart_fit_second")
heart <- adequacy(third, 2)
saveRDS(heart, file = "arbutus/heart_fit_third")

kidney <- adequacy(first, 3)
saveRDS(kidney, file = "arbutus/kidney_fit_first")
kidney <- adequacy(second, 3)
saveRDS(kidney, file = "arbutus/kidney_fit_second")
kidney <- adequacy(third, 3)
saveRDS(kidney, file = "arbutus/kidney_fit_third")

testis <- adequacy(first, 4)
saveRDS(testis, file = "arbutus/testis_fit_first")
testis <- adequacy(second, 4)
saveRDS(testis, file = "arbutus/testis_fit_second")
testis <- adequacy(third, 4)
saveRDS(testis, file = "arbutus/testis_fit_third")

cerebellum <- adequacy(first, 5)
saveRDS(cerebellum, file = "arbutus/cerebellum_fit_first")
cerebellum <- adequacy(second, 5)
saveRDS(cerebellum, file = "arbutus/cerebellum_fit_second")
cerebellum <- adequacy(third, 5)
saveRDS(cerebellum, file = "arbutus/cerebellum_fit_third")

liver <- adequacy(first, 6)
saveRDS(liver, file = "arbutus/liver_fit_first")
liver <- adequacy(second, 6)
saveRDS(liver, file = "arbutus/liver_fit_second")
liver <- adequacy(third, 6)
saveRDS(liver, file = "arbutus/liver_fit_third")
