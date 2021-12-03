#run Arbutus

library(arbutus)
library(geiger)
library(tidyverse)

#Load test data
brain <- c(readRDS("arbutus/brain_fit_first"), readRDS("arbutus/brain_fit_second"), readRDS("arbutus/brain_fit_third"))
cerebellum <- c(readRDS("arbutus/cerebellum_fit_first"), readRDS("arbutus/cerebellum_fit_second"), readRDS("arbutus/cerebellum_fit_third"))
heart <- c(readRDS("arbutus/heart_fit_first"), readRDS("arbutus/heart_fit_second"), readRDS("arbutus/heart_fit_third"))
kidney <- c(readRDS("arbutus/kidney_fit_first"), readRDS("arbutus/kidney_fit_second"), readRDS("arbutus/kidney_fit_third"))
liver <- c(readRDS("arbutus/liver_fit_first"), readRDS("arbutus/liver_fit_second"), readRDS("arbutus/liver_fit_third"))
testis <- c(readRDS("arbutus/testis_fit_first"), readRDS("arbutus/testis_fit_second"), readRDS("arbutus/testis_fit_third"))

all <- c(brain, cerebellum, heart, kidney, liver, testis)
rm(brain, cerebellum, heart, kidney, liver, testis)

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
    unit_tree$phy["edge.length"] <- map(unit_tree$phy["edge.length"], rm_neg) #fixing negative edge lengths
    print(unit_tree$phy["edge.length"])
    obs <- calculate_pic_stat(unit_tree)
    sim_tree <- simulate_char_unit(unit_tree)
    sim <- calculate_pic_stat(sim_tree)
    arby[[count]] <- compare_pic_stat(obs,sim)
    count = count + 1
  }
}
arby
}

test <- run_arb(all[1:3])

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

p_piv <- p_df %>% pivot_longer(cols = everything(), names_to = "tstat")

p_piv %>% ggplot(aes(value)) + geom_histogram() + facet_wrap(~tstat)

ggsave("arbutus/pvals_all.png")
