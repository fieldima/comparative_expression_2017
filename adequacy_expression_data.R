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
  result <- data.frame(
    m.sig = vector("double", length = length(tree)),
    c.var = vector("double", length = length(tree)),
    s.var = vector("double", length = length(tree)),
    s.asr = vector("double", length = length(tree)),
    s.hgt = vector("double", length = length(tree)),
    d.cdf = vector("double", length = length(tree))
    )
  count <- 1
  for (tr in tree){
    phylo <- tr@phylo
    dat <- tr@data %>% select(label, Averaged.RPKM.brain:Averaged.RPKM.liver) %>%
      slice(1:((nrow(tr@data)+1)/2)) %>% column_to_rownames(var = "label") %>%
      as.matrix()
    dat[is.na(dat)] = 0
    arb <- arbutus(phylo, data = dat[,part])$p.values
    result$m.sig <- arb[1]
    result$c.var <- arb[2]
    result$s.var <- arb[3]
    result$s.asr <- arb[4]
    result$s.hgt <- arb[5]
    result$d.cdf <- arb[6]
    count <- count + 1
  }
  result
}

#Measure adequacy
brain <- adequacy(calibrated, 1)
heart <- adequacy(calibrated, 2)
kidney <- adequacy(calibrated, 3)
testis <- adequacy(calibrated, 4)
cerebellum <- adequacy(calibrated, 5)
liver <- adequacy(calibrated, 6)

save.image(file = "arbutus/adequacy.RData")
