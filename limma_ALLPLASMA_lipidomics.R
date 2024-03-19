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
library(janitor)

setwd(dir = "chg306/INFLAMMATION_TRANSFER/INFLAMMATION_TRANSFER/")

##maternal plasma1
trans_df <- read_csv("../data/lipidomics/mat_plasma1_format.csv") %>%
  as.data.frame() 

#split exposure and timepoint and format
split_data <- strsplit(trans_df$exposure, "_")

# Extract 'exposure' and 'timepoint'
trans_df$exposure <- sapply(split_data, function(x) x[1])
trans_df$timepoint <- as.integer(sapply(split_data, function(x) x[2]))

#corrected
formatted_df <-
trans_df %>%
  relocate(timepoint, .after = exposure) %>%
  mutate(timepoint=factor(timepoint, levels=c("2", "5", "12", "24")))  %>%
  mutate(exposure=factor(exposure))

formatted_df$sample_id <-gsub("P","",formatted_df$sample_id)
    
#coldata for limmamodel
coldata <-
  formatted_df %>%
  column_to_rownames("sample_id") %>%
  .[,1:2]
  
## make tibble with sample id's as columns and make them numeric
data_df <- 
  formatted_df %>% .[-c(2:3)] %>%
  t() %>%
  row_to_names(row_number = 1)

# Extract row names to save after conversion to numeric
row_names <- rownames(data_df)

convert_number <- function(x){
  x <- as.character(x)
  #x <- gsub(pattern = ",", replacement = ".",x = x, fixed = TRUE)
  x <- as.numeric(x)
  return(x)
}

final_df <- 
  apply(data_df, 2, convert_number) 

rownames(final_df) <- row_names

##Go through each row and determine if a value is zero
row_sub = apply(final_df, 1, function(row) all(row !=0 ))
##Subset as usual
final_df <- final_df[row_sub,]


## log2 and scaling
data_norm <-
  scaling(final_df, type = "pareto")


#individual plots
## format for individual plots
data_long_values <-
  data_norm %>%
  as.data.frame() %>%
  rownames_to_column("metabolite") %>%
  pivot_longer(2:45,
               names_to = "sample_id",
               #names_repair = "check_unique",
               values_to = "peak_intensity",
  )

coldata_col <-
  coldata %>%
  rownames_to_column("sample_id")

joined_data <-
  inner_join(data_long_values, coldata_col, by = "sample_id")

joined_data %>%
  write_csv("data/lipidomics/USE_pareto_maternal_plasma_lipidomics_data.csv")



##go on with limma formatting
data_log <- 
  log2(data_norm + 1)

##design and run limma model

##combine timepoint and exposure
combo <- paste( coldata$exposure,coldata$timepoint,sep="_")

####
##make limma design
design <- model.matrix( ~ 0 + combo)
colnames(design) <- gsub("combo","", colnames(design)) 

## make contrasts
contrasts <- makeContrasts(
  W2 = lps_2 - ctl_2,
  W5 = lps_5 - ctl_5,
  W12 = lps_12 - ctl_12,
  W24 = lps_24 - ctl_24,
  levels = design)

##make limmafit and ebayes stats
fit <- lmFit(data_log, design)
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

## all limma, USE
all %>%
  write_csv("data/lipidomics/limma_maternal_plasma1_lipidomics.csv")

#with simplified names column too
all$metabolite_simple <- 
  sub("^(.*[0-9]+:[0-9]+).*", "\\1", all$metabolite)
  
all_simple <-
  all %>%
  relocate(metabolite_simple, .after = metabolite) %>%
  write_csv("data/lipidomics/limma_maternal_plasma1_lipidomics_simple.csv")

##significant FDR values
sign <-
all_simple %>%
  filter(adj.P.Val < 0.05) %>%
  distinct(metabolite)

## top50 fdr
all %>%
  head(n=50) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
write_csv("data/lipidomics/maternal_plasma1_top50_fdr.csv")

## top100 fdr
all_simple %>%
  head(n=100) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/maternal_plasma1_top100_fdr.csv")







### MATERNAL PLASMA 2
trans_df <- read_csv("data/lipidomics/mat_plasma2_format.csv")

trans_df[19, 2] <- 'ctl_12'
trans_df[20, 2] <- 'lps_9'
trans_df[21, 2] <- 'ctl_9'
trans_df[26, 2] <- 'lps_9'
trans_df[2, 2] <- 'ctl_9'
trans_df[7, 2] <- 'ctl_6'
trans_df[8, 2] <- 'lps_6'
trans_df[14, 2] <- 'ctl_6'
trans_df[15, 2] <- 'lps_6'
trans_df[16, 2] <- 'ctl_9'
trans_df[17, 2] <- 'lps_9'

trans_df %>%
  write_csv("data/lipidomics/mat_plasma2_format_corrected.csv")

#split exposure and timepoint and format
split_data <- strsplit(trans_df$exposure, "_")

# Extract 'exposure' and 'timepoint'
trans_df$exposure <- sapply(split_data, function(x) x[1])
trans_df$timepoint <- as.integer(sapply(split_data, function(x) x[2]))

#corrected
formatted_df <-
  trans_df %>%
  relocate(timepoint, .after = exposure) %>%
  mutate(timepoint=factor(timepoint, levels=c("6", "9", "12")))  %>%
  mutate(exposure=factor(exposure)) %>%
  mutate(sample_id = sub("M", "MP", sample_id))

#coldata for limmamodel
coldata <-
  formatted_df %>%
  column_to_rownames("sample_id") %>%
  .[,1:2]

## make tibble with sample id's as columns and make them numeric
data_df <- 
  formatted_df %>% .[-c(2:3)] %>%
  t() %>%
  row_to_names(row_number = 1)

# Extract row names to save after conversion to numeric
row_names <- rownames(data_df)

# Convert all columns to numeric
numeric_df <-
  as.data.frame(apply(data_df,    
                      2,
                      function(x) as.numeric(as.character(x))))

# Add row names as a new column in df2
numeric_df$RowNames <- row_names 

final_df <-
  numeric_df %>%
  column_to_rownames("RowNames")

##Go through each row and determine if a value is zero
row_sub = apply(final_df, 1, function(row) all(row !=0 ))
##Subset as usual
final_df <- final_df[row_sub,]


## log2 and scaling
data_norm <-
  scaling(final_df, type = "pareto")

data_log <- 
  log2(data_norm + 1)

##design and run limma model

##combine timepoint and exposure
combo <- paste( coldata$exposure,coldata$timepoint,sep="_")

####
##make limma design
design <- model.matrix( ~ 0 + combo)
colnames(design) <- gsub("combo","", colnames(design)) 

## make contrasts
contrasts <- makeContrasts(
  W6 = lps_6 - ctl_6,
  W9 = lps_9 - ctl_9,
  W12 = lps_12 - ctl_12,
  levels = design)

##make limmafit and ebayes stats
fit <- lmFit(data_log, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

results_w6 <- topTable(efit, coef=c("W6") ,Inf) %>%
  add_column(timepoint = 6) %>%
  rownames_to_column(var = "metabolite")
results_w9 <- topTable(efit, coef=c("W9") ,Inf) %>%
  add_column(timepoint = 9) %>%
  rownames_to_column(var = "metabolite")
results_w12 <- topTable(efit, coef=c("W12") ,Inf) %>%
  add_column(timepoint = 12) %>%
  rownames_to_column(var = "metabolite")

#make tibble of them all
all <-
  rbind.data.frame(results_w6, results_w9, results_w12) %>%
  #rownames_to_column(var = "metabolite") %>%
  .[order(.$adj.P.Val),] %>%
  mutate(timepoint=factor(timepoint, levels=c("6", "9", "12")))  

## all limma, USE
all %>%
  write_csv("data/lipidomics/limma_maternal_plasma2_lipidomics.csv")

##significant FDR values
sign <-
  all %>%
  filter(adj.P.Val < 0.05)

## top50 fdr
all %>%
  head(n=56) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/maternal_plasma2_top50_fdr.csv")

## top100 fdr
all %>%
  head(n=100) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/maternal_plasma2_top100_fdr.csv")





## FETAL PLASMA

trans_df <- read_csv("data/lipidomics/fet_plasma2_format.csv")

trans_df[19, 2] <- 'lps_9'
trans_df[20, 2] <- 'ctl_9'
trans_df[25, 2] <- 'lps_9'
trans_df[1, 2] <- 'ctl_9'
trans_df[6, 2] <- 'ctl_6'
trans_df[7, 2] <- 'lps_6'
trans_df[13, 2] <- 'ctl_6'
trans_df[14, 2] <- 'lps_6'
trans_df[15, 2] <- 'ctl_9'
trans_df[16, 2] <- 'lps_9'
trans_df[18, 2] <- 'lps_12'

trans_df %>%
  write_csv("data/lipidomics/fet_plasma2_format_corrected.csv")



#split exposure and timepoint and format
split_data <- strsplit(trans_df$exposure, "_")

# Extract 'exposure' and 'timepoint'
trans_df$exposure <- sapply(split_data, function(x) x[1])
trans_df$timepoint <- as.integer(sapply(split_data, function(x) x[2]))

#corrected
formatted_df <-
  trans_df %>%
  relocate(timepoint, .after = exposure) %>%
  mutate(timepoint=factor(timepoint, levels=c("6", "9", "12")))  %>%
  mutate(exposure=factor(exposure)) %>%
  mutate(sample_id = sub("F", "FP", sample_id))

#coldata for limmamodel
coldata <-
  formatted_df %>%
  column_to_rownames("sample_id") %>%
  .[,1:2]

## make tibble with sample id's as columns and make them numeric
data_df <- 
  formatted_df %>% .[-c(2:3)] %>%
  t() %>%
  row_to_names(row_number = 1)

# Extract row names to save after conversion to numeric
row_names <- rownames(data_df)

# Convert all columns to numeric
numeric_df <-
  as.data.frame(apply(data_df,    
                      2,
                      function(x) as.numeric(as.character(x))))

# Add row names as a new column in df2
numeric_df$RowNames <- row_names 

final_df <-
  numeric_df %>%
  column_to_rownames("RowNames")

##Go through each row and determine if a value is zero
row_sub = apply(final_df, 1, function(row) all(row !=0 ))
##Subset as usual
final_df <- final_df[row_sub,]


## log2 and scaling
data_norm <-
  scaling(final_df, type = "pareto")

data_log <- 
  log2(data_norm + 1)

##design and run limma model

##combine timepoint and exposure
combo <- paste( coldata$exposure,coldata$timepoint,sep="_")

####
##make limma design
design <- model.matrix( ~ 0 + combo)
colnames(design) <- gsub("combo","", colnames(design)) 

## make contrasts
contrasts <- makeContrasts(
  W6 = lps_6 - ctl_6,
  W9 = lps_9 - ctl_9,
  W12 = lps_12 - ctl_12,
  levels = design)

##make limmafit and ebayes stats
fit <- lmFit(data_log, design)
fit2 <- contrasts.fit(fit, contrasts=contrasts)
efit <- eBayes(fit2)

results_w6 <- topTable(efit, coef=c("W6") ,Inf) %>%
  add_column(timepoint = 6) %>%
  rownames_to_column(var = "metabolite")
results_w9 <- topTable(efit, coef=c("W9") ,Inf) %>%
  add_column(timepoint = 9) %>%
  rownames_to_column(var = "metabolite")
results_w12 <- topTable(efit, coef=c("W12") ,Inf) %>%
  add_column(timepoint = 12) %>%
  rownames_to_column(var = "metabolite")

#make tibble of them all
all <-
  rbind.data.frame(results_w6, results_w9, results_w12) %>%
  #rownames_to_column(var = "metabolite") %>%
  .[order(.$adj.P.Val),] %>%
  mutate(timepoint=factor(timepoint, levels=c("6", "9", "12")))  

## all limma, USE
all %>%
  write_csv("data/lipidomics/limma_fetal_plasma2_lipidomics.csv")

##significant FDR values
sign <-
  all %>%
  filter(adj.P.Val < 0.05)

## top50 fdr
all %>%
  head(n=56) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/fetal_plasma2_top50_fdr.csv")

## top100 fdr
all %>%
  head(n=100) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/fetal_plasma2_top100_fdr.csv")

  
# normality - skip for now
  ggdensity(joined_data$peak_intensity, xlim = c(0.0,6))
  
  shapiro.test(joined_data$peak_intensity)
  
  ks.test(joined_data$peak_intensity,
          alternative = c("two.sided"),
          exact = NULL, simulate.p.value = FALSE, B = 2000)

  
  
  
