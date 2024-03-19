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

lipid_df <- read_csv("../data/lipidomics/use_maternal_liver_raw.csv")

## format df with sample names, timepoints, exposure
lipid_df$Name[1] <- 'exposure'

trans_df <- lipid_df %>%
    as.data.frame() %>%
    column_to_rownames("Name") %>%
    t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") 

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

##get no of replicates for plots, make table

repl_no <- 
  formatted_df %>%
  select(sample_id, exposure, timepoint) %>%
  group_by(exposure, timepoint) %>%
  dplyr::summarise(replicates = n()) 


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
numeric_df <- ##changed name from numeric_df
as.data.frame(apply(data_df,    
                    2,
                    function(x) as.numeric(as.character(x))))

# add tissue id's sto sample names (SKIP if creating INDIVIDUAL LIPIDS csv file)
new_names <- paste("MLi", names(numeric_df), sep = "")
names(numeric_df) <- new_names

# Add row names as a new column in df2 - to keep metabolite 
numeric_df$RowNames <- row_names 

final_df <-
  numeric_df %>%
  column_to_rownames("RowNames")

##Go through each row and determine if a value is zero
row_sub = apply(final_df, 1, function(row) all(row !=0 ))
##Subset as usual
final_df <- final_df[row_sub,]

#dont run lines below
### write data for Mixomics
write.table(rbind.data.frame(trans_df$exposure,final_df),
            "../data/lipidomics/use_maternal_liver_raw_formatted.csv", quote=FALSE, row.names=FALSE, sep=",")

## format for individual plots
data_long_values <-
data_norm %>%
  rownames_to_column("metabolite") %>%
  pivot_longer(2:57,
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
  write_csv("data/lipidomics/USE_pareto_maternal_liver_lipidomics_data.csv")


##end of lines not to run, continue for limma

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
  W2 = lps_2 - ctr_2,
  W5 = lps_5 - ctr_5,
  W12 = lps_12 - ctr_12,
  W24 = lps_24 - ctr_24,
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
  write_csv("data/lipidomics/limma_maternal_liver_lipidomics.csv")

#with simplified names column too
all$metabolite_simple <- 
  sub("^(.*[0-9]+:[0-9]+).*", "\\1", all$metabolite)

all_simple <-
  all %>%
  relocate(metabolite_simple, .after = metabolite) %>%
  write_csv("data/lipidomics/limma_maternal_liver_lipidomics_simple.csv")

##significant FDR values
sign <-
all %>%
  filter(adj.P.Val < 0.05)

## top50 fdr
all %>%
  head(n=51) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
write_csv("data/lipidomics/maternal_liver_top50_fdr.csv")

## top100 fdr
all_simple %>%
  head(n=100) %>%
  dplyr::pull(metabolite) %>%
  tibble() %>%
  write_csv("data/lipidomics/maternal_liver_top100_fdr.csv")

  
# normality - skip for now
  ggdensity(joined_data$peak_intensity, xlim = c(0.0,6))
  
  shapiro.test(joined_data$peak_intensity)
  
  ks.test(joined_data$peak_intensity,
          alternative = c("two.sided"),
          exact = NULL, simulate.p.value = FALSE, B = 2000)
