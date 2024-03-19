.libPaths(c("/maps/projects/scarball/people/chg306/R_libs",
            "/opt/software/R/4.2.2/lib64/R/library",
            "/maps/projects/scarball/apps/Modules/software/R-packages/4.2.2"))

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)
library(Cairo)
library(MetabolAnalyze)


lipid_df <- read_csv("data/lipidomics/use_data_raw.csv")

## create tibble with sample names, timepoints, exposure
coldata <- lipid_df %>%
    as.data.frame() %>%
    .[1:2,] %>%
    remove_rownames() %>% 
    column_to_rownames("metabolite") %>% 
    t %>%
    as.data.frame %>% 
    mutate(timepoint=factor(timepoint, levels=c("2", "5", "12", "24")))  %>%
    mutate(exposure=factor(exposure)) 
    

## make tibble with sample id's as columns and make them numeric
data_df <- lipid_df %>% as.data.frame %>% .[-c(1:2),] %>%
    remove_rownames() %>% 
    column_to_rownames("metabolite")
for (i in c(1:ncol(data_df)))
{
   data_df[, i] <-  as.numeric(data_df[,i])
}

##Go through each row and determine if a value is zero
row_sub = apply(data_df, 1, function(row) all(row !=0 ))
##Subset as usual
data_df <- data_df[row_sub,]

##combine timepoint and exposure
combo <- paste( coldata$exposure,coldata$timepoint,sep="_")

##make limma design
design <- model.matrix( ~ 0 + combo)
colnames(design) <- gsub("combo","", colnames(design)) 

## make contrasts
contrasts <- makeContrasts(
    W2 = lps_2 - ctr_2,
    W5 = lps_5 - ctr_5,
    W12 = lps_12 - ctr_12,
    W24 = lps_24 - ctr_24,
    levels = design)

##first scale data (also use for individual lipid plots, but then don't transform them afterwards)
data_pareto <-
scaling(data_df, type = "pareto")

## use log2 to transform for limma, OMVENDT
data_transformed <- 
  log2(data_pareto + 1)

#format data further for limma
data_long_values <-
  data_transformed %>%
rownames_to_column("metabolite") %>%
  pivot_longer(2:77,
    names_to = "sample_id",
    #names_repair = "check_unique",
    values_to = "peak_intensity",
  )

coldata_col <-
coldata %>%
  rownames_to_column("sample_id")

joined_data <-
inner_join(data_long_values, coldata_col, by = "sample_id") 

# normality
ggdensity(joined_data$peak_intensity, xlim = c(0.0,6))

shapiro.test(joined_data$peak_intensity)

ks.test(joined_data$peak_intensity,
        alternative = c("two.sided"),
        exact = NULL, simulate.p.value = FALSE, B = 2000)


## for individual plots
joined_data %>%
  write_csv("data/lipidomics/USE_pareto_fetal_liver_lipidomics_data.csv")

##make limmafit and ebayes stats
fit <- lmFit(data_transformed, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

results_w2 <- topTable(efit, coef=c("W2") ,Inf) %>%
  add_column(timepoint = 2) %>%
rownames_to_column(var = "metabolite")
results_w5 <- topTable(efit, coef=c("W5") ,Inf) %>%
  add_column(timepoint = 5) %>%
  rownames_to_column(var = "metabolite")
results_w12 <- topTable(efit, coef=c("W12") ,Inf) %>%
  add_column(timepoint = 12) %>%
  rownames_to_column(var = "metabolite")
results_w24 <- topTable(efit, coef=c("W24") ,Inf) %>%
  add_column(timepoint = 24) %>%
  rownames_to_column(var = "metabolite")

#make tibble of them all
all <-
  rbind.data.frame(results_w2, results_w5, results_w12, results_w24) %>%
  #rownames_to_column(var = "metabolite") %>%
  .[order(.$adj.P.Val),] %>%
  mutate(timepoint=factor(timepoint, levels=c("2", "5", "12", "24")))  

## all limma, USE - or use version with necessary columns further down
all%>%
  write_csv("data/lipidomics/limma_lipidomics.csv")

#scaled
?

##significant FDR values
sign <-
all %>%
  filter(adj.P.Val < 0.05)

## top50 fdr
all %>%
  head(n=56) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
write_csv("data/lipidomics/top50_fdr.csv")

## top100 fdr
all %>%
  head(n=100) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/top100_fdr.csv")


## Only have useful columns
t <-
all %>%
 # left_join(all, lps_tibble, by = c("timepoint", "metabolite")) %>%
#.[,c(1:3,6,8:11)] %>%
  .[,c(1:2,6,8)]  %>%
write_csv("data/lipidomics/use_limma_fetal_liver.csv")
