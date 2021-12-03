#Diagnose

#Load packages
library(ape)
library(geiger)
library(arbutus)
library(treeio)
library(ggtree)
library(tidyverse)

#Load treedata
calibrated <- readRDS("export/calibrated_trees") %>% discard(~ Nnode2(.x) < 9)
test <- calibrated[1:1000]

rm_neg <- function (edge){
  ifelse(edge < 0, 1e-06, edge)
}
cal <- vector("list", length = length(calibrated))
ult <- vector(mode = "logical", length = length(calibrated))
count = 1
for(cali in calibrated){
  #print(cali@data)
  #print(names(cali@data))
  #print(cali@phylo$tip.label)
  dat <- cali@data %>% filter(!is.na(Mean))
  dv <- dat$Mean
  names(dv) <- dat$label
  tr <- geiger::treedata(cali@phylo, dv)
  unit.tree <- make_unit_tree(tr$phy, data = tr$data[,1])
  unit.tree$phy["edge.length"] <- map(unit.tree$phy["edge.length"], rm_neg)
  ult[count] <- is.ultrametric(unit.tree$phy)
  count = count +1
  #obs <- calculate_pic_stat(unit.tree)
  #print(obs)
  #sim.dat <- simulate_char_unit(unit.tree)
  #print(sim.dat)
  #sim <- calculate_pic_stat(sim.dat)
  #print(sim)
  #res <- compare_pic_stat(obs, sim)
  #print(res)
}


#Removal of negatives caused 6 trees to no longer be ultrametric
neg <- calibrated[!ult]
pos <- calibrated[ult]
