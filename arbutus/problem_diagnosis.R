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
test <- calibrated[1:10]

rm_neg <- function (edge){
  ifelse(edge < 0, 1e-06, edge)
}
result <- vector("list", length = 6*length(calibrated))
count = 1
for(cali in calibrated){
  dat <- cali@data %>% filter(!is.na(Mean)) %>% select(label,Averaged.RPKM.brain:Averaged.RPKM.liver, Mean)
  cali@phylo["edge.length"] <- map(cali@phylo["edge.length"], rm_neg)
  dat[is.na(dat)] = 0.00001
  dbrain <- dat[,2]  
  names(dbrain) <- dat$label
  dheart <- dat[,3]  
  names(dheart) <- dat$label
  dkidney <- dat[,4]  
  names(dkidney) <- dat$label
  dtestis <- dat[,5]  
  names(dtestis) <- dat$label
  dcerebellum <- dat[,6]  
  names(dcerebellum) <- dat$label
  dliver <- dat[,7]  
  names(dliver) <- dat$label
  dl <- list(dbrain, dheart, dkidney, dtestis, dcerebellum, dliver)
  dv <- dat[,8]
  names(dv) <- dat$label
  g.tr <- geiger::treedata(cali@phylo, dv)
  g.tree <- make_unit_tree(g.tr$phy, data = g.tr$data[,1])
  g.sim <- simulate_char_unit(g.tree) %>% calculate_pic_stat()
  for(d in dl){
    tr <- geiger::treedata(cali@phylo, d)
    unit.tree <- make_unit_tree(tr$phy, data = tr$data[,1])
    obs <- calculate_pic_stat(unit.tree)
    res <- compare_pic_stat(obs, g.sim)
    result[[count]] = res
    count = count + 1
  }
}

pvals <- map(result, ~ .x$p.values)

arbutus_transform <- function ( pval , tib) {
  len <- length(pval)
  df <- t(data.frame(pval))
  if(missing(tib)) tib = FALSE
  ifelse(tib == TRUE, {
    fin <- as_tibble(df)
  }
  , fin <- data.frame(df))
  row.names(fin) <- c(1:len)
  fin
}
p_df <- arbutus_transform(pvals)

p_piv <- p_df %>% pivot_longer(cols = everything(), names_to = "tstat")

saveRDS(p_piv, file = "arbutus/samesimpiv")
p_piv %>% ggplot(aes(value)) + geom_histogram() + facet_wrap(~tstat)
ggsave("arbutus/samesimpvals.png")

#Removal of negatives caused 6 trees to no longer be ultrametric
#neg <- calibrated[!ult]
#pos <- calibrated[ult]