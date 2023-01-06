#One big script with everything

#fitContinuous on comparative expression data

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


#Load data
brain <- c(readRDS("arbutus/brain_fit_first"), readRDS("arbutus/brain_fit_second"), readRDS("arbutus/brain_fit_third"))
cerebellum <- c(readRDS("arbutus/cerebellum_fit_first"), readRDS("arbutus/cerebellum_fit_second"), readRDS("arbutus/cerebellum_fit_third"))
heart <- c(readRDS("arbutus/heart_fit_first"), readRDS("arbutus/heart_fit_second"), readRDS("arbutus/heart_fit_third"))
kidney <- c(readRDS("arbutus/kidney_fit_first"), readRDS("arbutus/kidney_fit_second"), readRDS("arbutus/kidney_fit_third"))
liver <- c(readRDS("arbutus/liver_fit_first"), readRDS("arbutus/liver_fit_second"), readRDS("arbutus/liver_fit_third"))
testis <- c(readRDS("arbutus/testis_fit_first"), readRDS("arbutus/testis_fit_second"), readRDS("arbutus/testis_fit_third"))

all <- c(brain, cerebellum, heart, kidney, liver, testis)
rm(brain, cerebellum, heart, kidney, liver, testis)

#First need to pull out model names
model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- unlist(f)
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  data.frame(OU = ou, BM = bm, EB = eb)
}

df <- model_count(all)
saveRDS(all, file = "arbutus/fitdata")

b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col()
ggsave(filename = "arbutus/AIC_plot.png")


rm_neg <- function (edge){
  ifelse(edge < 0, 1e-06, edge)
}

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(fit in fits){
    for(f in fit){
      class(f) <- "gfit"
      unit_tree <- make_unit_tree(f)
      unit_tree$phy["edge.length"] <- purrr::map(unit_tree$phy["edge.length"], rm_neg) #fixing negative edge lengths
      print(paste(count, "done"))
      obs <- calculate_pic_stat(unit_tree)
      sim_tree <- simulate_char_unit(unit_tree)
      sim <- calculate_pic_stat(sim_tree)
      arby[[count]] <- compare_pic_stat(obs,sim)
      count = count + 1
    }
  }
  arby
}

arb_result <- run_arb(all)
saveRDS(arb_result, file = "arbutus/arbutus_results_all")

count = 1
pvals <- vector("list", length = length(arb_result))
for(arb in arb_result){
  pvals[[count]] = arb$p.values
  count = count + 1
}

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
saveRDS(object = p_df, file = "arbutus/p_vals_df")

#remove m.sig
p_piv <- p_df %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat")

#reorder facets
p_piv$tstat <- factor(p_piv$tstat, levels = c("c.var", "s.var", "s.asr", "s.hgt", "d.cdf"))
p_piv %>% ggplot(aes(value, after_stat(density))) + geom_histogram() + facet_wrap(~tstat, nrow = 1) + theme_bw()

ggsave("arbutus/pvals_all.png")