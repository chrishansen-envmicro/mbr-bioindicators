#load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(phyloseq)
library(wesanderson)
library(vegan)
library(indicspecies)

# contents: methods to perform, modify, visualize the output of Sperlea et al.'s (2021)
# supervised machine learning framework. This pipeline was applied to assess
# the bioindicator potential of bacterial communities on glacial sourcing (fi),
# relative to other physical and geochemical parameters.

## references

# Miller, J. B., Frisbee, M. D., Hamilton, T. L., & Murugapiran, S. K. (2021). 
# Recharge from glacial meltwater is critical for alpine springs and their microbiomes. 
# Environmental Research Letters, 16(6), 64012-. https://doi.org/10.1088/1748-9326/abf06b 

# Sperlea, T., Kreuder, N., Beisser, D., Hattab, G., Boenigk, J., & Heider, D. (2021). 
# Quantification of the covariation of lake microbiomes and environmental variables 
# using a machine learning‐based framework. Molecular Ecology, 30(9), 2131–2144. 
# https://doi.org/10.1111/mec.15872

#### Phyloseq Object Setup ####

#Import data from MSI: /panfs/roc/groups/6/hamil689/shared/sen/JM_second_paper/res/16SII
sharedfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated_4.12.csv"


#Returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

#Import metadata -> remove blanks and samples lacking recharge_source
map <- read.csv(mapfile)
# rename location values to improve figure clarity
map$location <- recode(map$location, MH = "Mount Hood") %>% 
  recode(GNP = "Glacier National Park") 

map.recharge <- filter(map, ncbi_name_18s != "") # remove samples without 18S data
map.recharge <- sample_data(map.recharge) # format sample data
rownames(map.recharge) <- map.recharge$ncbi_name_18s # align rownames across metadata .tax and .shared files

# merge morthur files (.taxonomy, .shared) with metadata
moth_merge <- merge_phyloseq(mothur_data, map.recharge)

#replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# remove springs and streams
stream.spring <- c("Stream", "Spring")
data_spst <- subset_samples(moth_merge, sample_type %in% stream.spring)


#remove singletons
data_filt <- filter_taxa(data_spst, function (x) {sum(x > 0) > 1}, prune=TRUE)

#develop rarefaction curves and visualize distribution of read lengths: OPTIONAL
tab <- otu_table(data_filt)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
library(vegan)
rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)
abline(v=14000)

# how many OTUs in each sample?
df <- tab %>% 
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")
otu_summary <- df %>%
  group_by(SampleID) %>%
  summarize(non_zero_otus = sum(across(where(is.numeric), ~ .x != 0)))

# how many reads in each sample?
row_sums <- rowSums(tab)
sorted_row_sums <- sort(row_sums)
print(sorted_row_sums)

#identify which samples had fewer than 10^3 reads and remove (12)
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (20000))
min(sample_sums(data_filt_clean)) #check min at 11000. revise above to be 90% of this min (11919)
sort(sample_sums(data_filt))

#how many samples remain from each mountain range?
sample_data(data_filt_clean) %>%
  group_by(location) %>%
  summarise(n = n())

#hist(sample_sums(data_spst), main="Histogram: Read Counts", xlab="Total Reads", 
#border="blue", col="green", las=1, breaks=12)

#Rarefy the samples without replacement. Rarefaction is used to simulate even number 
#of reads per sample. In this example, the rarefaction depth chosen is the 90% of the 
#minimum sample depth in the dataset (in this case 459 reads per sample).
#rarefy based on lowers read count about 1k. (1015)
#https://micca.readthedocs.io/en/latest/phyloseq.html
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=min(sample_sums(data_filt_clean)), 
                                       replace=F)
#perform srs to standardize counts: https://rdrr.io/github/jfq3/QsRutils/man/srs_p.html
#data_filt_clean_srs <- srs_p(data_filt_clean)

sum(otu_table(data_filt_rarified)) # n_reads
nsamples(data_filt_rarified) # n_samples
ntaxa(data_filt_rarified) # n_OTUs   

# split MH and GNP samples
data.filt.rarified.mh <- subset_samples(data_filt_rarified, location == "Mount Hood")
data.filt.rarified.gnp <- subset_samples(data_filt_rarified, location == "Glacier National Park")

#### Sperlea Covariation Template (https://onlinelibrary.wiley.com/doi/10.1111/mec.15872) ####

# This is sample code for the covariation framework as described in 
#  Sperlea et al. (2021), "Quantification of the covariation of lake microbiomes and
#  environmental variables using a machine learning-based framework", Molecular Ecology
#  Needs to be adjusted for the datasets used.

library("data.table")
library(caret)
library(indicspecies)

args = commandArgs(trailingOnly=TRUE)

##Arguments:
# 1 - environmental data
# 2 - taxa / OTU table
# 3 - path to the folder to save models to (as RDS files)
# 4 - path to the folder to save model predictions to (as csv files)
# 5 - path to the folder to save model evaluations to (as csv files)
# 6 - target feature (as string)
# 7 - machine learning model (as string that will be passed on to the train() 
#     function from caret as "method" parameter

###### START Function Creation ######

fcbf_filter <- function(in_data, target_data, method, threshold){
  ##Fast correlaton-based filter as described by Yu & Liu (2003)
  tmp.data = in_data 
  features = c()
  
  vars = apply(tmp.data, 2, var)
  dels = which(vars == 0)
  
  if(length(dels) > 0){
    tmp.data = tmp.data[,-dels]
  }
  tmp.data.names = colnames(tmp.data)
  
  counter = 0
  while(!is.null(dim(tmp.data))){
    cor.class = abs(as.vector(cor(tmp.data, target_data, method=method)))
    imp = t(as.matrix(cor.class))
    best = tmp.data.names[which.max(as.vector(imp[1,]))]
    features = c(features, best)
    
    ##remove everything that correlates with best above cutoff
    cor.best = abs(as.vector(cor(tmp.data, tmp.data[, best], method=method)))
    w = which(cor.best > threshold)
    
    tmp.data = tmp.data[,-w]
    tmp.data.names = tmp.data.names[-w]
    counter = counter + 1
  }
  
  return(features)
}

#### Determine tertiles of variable in question ####

#### Complete this and update indval_filter accordingly after getting to final run stage ####
quantile(environmental_data$tds, c(0:3/3))

indval_filter <- function(in_data, target_data, cutoff = 0.05){
  ##Feature selection based on bioindicator analysis as described by 
  # Dufrene & Legendre, 1997. The target variable is split into tertiles by value
  # and the taxon / OTU table is preprocessed using Hellinger transformation
  
  dummy = target_data
  low_pos = which(dummy < 39) #lower tertile. MUST update for every env variable
  mid_pos = which(dummy >= 39) #middle tertile. MUST update for every env variable
  high_pos = which(dummy > 76) #higher tertile. MUST update for every env variable
  
  dummy[low_pos] = "low" 
  dummy[mid_pos] = "mid" 
  dummy[high_pos] = "high" 
  
  hellinger_data = as.data.frame(t(apply(in_data, 1, function(x) sqrt(x / sum(x)))))
  indval_output = multipatt(hellinger_data, dummy, control = how(nperm=999), func = "IndVal")$sign
  selected_variables = indval_output[indval_output$p.value < cutoff & indval_output$index != 7 & !(is.na(indval_output$p.value)),  , drop=F]
  
  return(rownames(selected_variables))
}

train_covariation_framework <- function(in_dataset, out_dataset, target_feature, 
                                        model, model_target, prediction_target, evaluation_target, 
                                        num_folds = 10, feature_selection_method){
  
  ##Create folds for cross-validation
  folds = createFolds(out_dataset[[target_feature]], num_folds)  
  
  collected_pred = c()
  collected_true = c()
  
  foldcounter = 1
  for (thisfold in folds){
    print(foldcounter)
    
    #Feature selection using either FCBF or indVal 
    if (feature_selection_method == "FCBF"){
      selected_variables = fcbf_filter(in_dataset[-thisfold,], 
                                       out_dataset[-thisfold, target_feature], "pearson", 0.6)
    }else if(feature_selection_method == "indVal"){
      selected_variables = indval_filter(in_dataset[-thisfold,], 
                                         out_dataset[-thisfold, target_feature])
    }
    
    ##create data subsets for training/prediction/evaluation
    train_in = in_dataset[-thisfold, selected_variables]
    train_out = out_dataset[-thisfold, target_feature]
    test_in = in_dataset[thisfold, selected_variables]
    test_out = out_dataset[thisfold, target_feature]
    
    ##train model
    print(model)
    if (model == "rf"){
      fit = train(train_in, train_out, method = model, importance = TRUE)
    }else{
      fit = train(train_in, train_out, method = model)
    }
    ##Save model for, e.g., feature importance analysis
    saveRDS(fit, file = paste(model_target, model, "_", target_feature, "_", 
                              as.character(foldcounter), ".RDS", sep = ""))
    
    ##predict environmental variable values based on the model and write to file
    pred <- predict(fit, test_in)
    pred_true = cbind(pred, test_out)
    colnames(pred_true) = c("pred", "true")
    write.csv(pred_true, paste(prediction_target, model, "_", target_feature, 
                               "_", as.character(foldcounter), ".csv", sep = ""))
    
    ##collect the predictions and the measured values for evaluation
    collected_pred = c(collected_pred, pred)
    collected_true = c(collected_true, test_out)
    
    foldcounter = foldcounter + 1
  }
  
  ##Calculate performance metrics and save to file
  evaluation = postResample(collected_pred, collected_true)
  write.csv(evaluation, paste(evaluation_target, model, "_", target_feature, 
                              ".csv", sep=""))
  
}

###### END Function Creation ######

# Both the environmental data as well as the taxonomic data/OTU tables need to 
# contain their samples in the same order and in rows; columns need to be named 
# and contain the different variables/taxa/OTUs. Neither can contain any NA 
# values or all-zero columns.

### Format OTU data (phyloseq -> indicspecies; https://github.com/joey711/phyloseq/issues/1212) ###
# Write out your phyloseq OTU table and export it

#STOP! Make sure you are using the correct rarified dataset (MH, GNP, Both)
write.csv(data.filt.rarified.gnp@otu_table,'otus_gnp.csv')
# Import phyloseq OTU table as an OTU table/dataframe
otu <-read.csv('otus_gnp.csv')
# do some shuffling of the OTU table
otu_flip <- as.data.frame(t(otu)) # makes it a dataframe and puts x into y and y into x (flips it)
names(otu_flip) <- as.matrix(otu_flip[1, ]) # renames columns
otu_flip<- otu_flip[-1, ] # removes first row
otu_flip_num<-as.data.frame(lapply(otu_flip, as.numeric)) # convert from character to number
otu_flip_num$ncbi_name_18s<-row.names(otu_flip) # places row names as sample ID column
# OK - now we have the OTU table fromatted correctly

#Repeat for metadata

#STOP! Make sure you are using the correct rarified dataset (MH, GNP, Both)
sam_dat <- as.matrix(sample_data(data.filt.rarified.gnp))
write.csv(sam_dat, "sam_dat.csv")
meta <- read.csv('sam_dat.csv')
## Join based on sample_id
otu_final<-left_join(meta, otu_flip_num , by = c("ncbi_name_18s" = "ncbi_name_18s")) # join based on sample IDs, assuming they're the same for both OTU table and metadata
#otu_final <- filter(otu_final, tds != "NA")
tax_data <- otu_final[,70:7880] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
hellinger_tax_data <- as.data.frame(t(apply(tax_data, 1, function(x) sqrt(x / sum(x)))))

# select the environmental parameters of interest
environmental_data <- select(otu_final, fi_best, elevation_m, ysi_temp, ysi_do, ysi_spc,
                             ysi_ph, ysi_orp, alkalinity_caco3, bicarbonate_hco3,
                             calcium, chloride, magnesium, potassium, sio2, sodium,
                             strontium, sulfate, tds, dist_from_summit)

### REP 1
# Specify output location and name for model, prediction, and evaluation files
modelfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/mod/mod_1"
predfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/pred/pred_1"
evalfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/eval_1"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "tds", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 10,
                            feature_selection_method = "indVal")
### REP 2
# Specify output location and name for model, prediction, and evaluation files
modelfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/mod/mod_2"
predfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/pred/pred_2"
evalfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/eval_2"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "tds", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 10,
                            feature_selection_method = "indVal")
### REP 3
# Specify output location and name for model, prediction, and evaluation files
modelfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/mod/mod_3"
predfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/pred/pred_3"
evalfolder <- "~/Desktop/Alpine_Hydro_Project/alpine_hydro/data/covar_evals_18S_6.9/gnp/eval_3"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "tds", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 10,
                            feature_selection_method = "indVal")

# use the indicspecies package to identify taxa that significant correlate with environmental variables
environmental_tert <- environmental_data %>%
  mutate(recharge_level = cut(fi_best,
                              quantile(fi_best, c(0:3/3)),
                              labels = c("low", "medium", "high"),
                              include.lowest = TRUE))
recharge_data <- environmental_tert$recharge_level
multi_test <- multipatt(x=hellinger_tax_data, cluster=recharge_data, control = how(nperm=999))
summary(multi_test)
#correct p-values?: https://stackoverflow.com/questions/55271042/r-indicspecies-package-multipatt-function-extract-values-from-summary-multip

# replace OTUs with species name
multipatt_output <- multi_test$sign
mp_filter <- filter(multipatt_output, p.value <= 0.05)

# extract tax data as a dataframe
write.csv(filter.data@tax_table,'tax.csv')

# import phyloseq OTU table as an OTU table/dataframe
tax <-read.csv('tax.csv')
mp_filter$X <- rownames(mp_filter) 
join <- left_join(mp_filter, tax)
join_filter <- filter(join, s.high == 1 & s.low == 0 & s.medium == 0)

# summarize number of bioindicators at various taxaonomic levels
join_filter %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  #top_n(n=7,wt = n) %>%
  arrange(Genus, desc(n))

### visualize covariance data --> "cov_data.csv" developed manually from evaluation outputs
cov_data <- read.csv("cov_data_18S_6.9.csv") %>%
  filter(Location == "Mount Hood" | Location == "Glacier National Park")

cov_mh <- filter(cov_data, Location == "Mount Hood") # extract MH data
cov_gnp <- filter(cov_data, Location == "Glacier National Park") # extract GNP data

cov_data_order <- cov_data %>% # reorder variables to reflect an increasing mean r_squared
  mutate(Env_Variable = fct_reorder(Env_Variable, Rsquared, mean, na.rm = TRUE,
                                    .desc = TRUE))
# plot
ggplot(cov_data_order, aes(x = Env_Variable, y = Rsquared, color = Location, shape = Location)) +
  geom_jitter(alpha = 0.8, width = 0.1, size = 5, show.legend = FALSE) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  geom_hline(yintercept = summarise(cov_mh, 
                                    mean(Rsquared, na.rm = TRUE)) %>% pull(), 
             lty = 2, color = wes_palette("Moonrise2")[2], lwd = 2)+ 
  geom_hline(yintercept = summarise(cov_gnp, 
                                    mean(Rsquared, na.rm = TRUE)) %>% pull(), 
             lty = 2, color = wes_palette("Moonrise2")[1], lwd = 2)+
  labs(y = bquote(~R^2)) +
  stat_summary(fun.data = "mean_cl_boot", #Add mean and 05% confidence intervals
               show.legend = TRUE,
               position = position_nudge(x = .3, y = 0),
               fun.args=(conf.int=0.95),
               size = 1,
               pch = 1) +
  coord_flip() +
  scale_y_continuous(limits = c(0,1)) +
  theme_light() +
  guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) +
  theme(axis.title = element_text(size = 20), #Modify text sizes to improve figure readability
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")

# statistical analysis for R2
t.test(Rsquared ~ Location, data = cov_data_order)
cov_sum <- cov_data_order %>%
  group_by(env_variable, cat) %>%
  summarise(mean = mean(r_squared))

# range of differences
cov_sum %>%
  group_by(env_variable) %>%
  summarise(diff = diff(mean)) %>%
  filter(env_variable == "Cl" |
           env_variable == "Na" |
           env_variable == "Sr" |
           env_variable == "K" |
           env_variable == "SiO2" |
           env_variable == "SO4" |
           env_variable == "DO" |
           env_variable == "Mg") %>%
  arrange(diff)
