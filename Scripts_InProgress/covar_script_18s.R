library(tidyverse)
library(readxl)
library(dplyr)
library(phyloseq)
library(wesanderson)
library(vegan)

#### Phyloseq Object Setup) ####

#Import data from MSI: /panfs/roc/groups/6/hamil689/shared/sen/JM_second_paper/res/16SII
sharedfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated.csv"

#Returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

#Import metadata -> remove blanks and samples lacking recharge_source
map <- read.csv(mapfile)
map$location <- recode(map$location, MH = "Mount Hood") %>% #rename location values to improve figure clarity
  recode(GNP = "Glacier National Park")

#filter to extract variables of interest and avoid duplicate env values from replicates
map_filt <- filter(map, sample_type == "Spring" | sample_type == "Stream") %>%
  filter(bio_rep != "NA") %>%
  filter(ncbi_name_18s != "")


map_filt <- sample_data(map_filt)
rownames(map_filt) <- map_filt$ncbi_name_18s


#Merge mothur data with metadata to create phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map_filt)

#Replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Remove non-bacterial mitochondria and chloropast seq. should be removed during amp processing (mother), but this is a double-check
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

#filter data (this is a choice depending on the project).
#below only takes otus with 1+ occurrences and <10% of total reads
filter.data <- filter_taxa(data, function (x) {sum(x > 0) > 1}, prune=TRUE)
#develop rarefaction curves and visualize distribution of read lengths: OPTIONAL
rarecurve(t(otu_table(filter.data)), step = 100, label = F)
abline(v=12000)

#identify which samples had fewer than 10^3 reads and remove (12)
data.filt.clean <- subset_samples(filter.data, sample_sums(filter.data) >= (6000))
min(sample_sums(data.filt.clean)) #check min at 8500. revise above to be 90% of this min (8792)

#how many samples remain from each mountain range?
sample_data(data.filt.clean) %>%
  group_by(location) %>%
  summarise(n = n())

#complete rarefaction
data.filt.rarified = rarefy_even_depth(data.filt.clean, rngseed=1, 
                                       sample.size=0.9*min(sample_sums(data.filt.clean)), 
                                       replace=F) %>%
  subset_samples(ysi_ph != "NA") #remove NAs for each subsample when attempting to address them
                 #calcium, chloride, magnesium, potassium, sio2, sodium,
                 #strontium, sulfate, tds)
  
#create equaivalents for MH and GNP
data.filt.rarified.mh <- subset_samples(data.filt.rarified, location == "Mount Hood")
data.filt.rarified.gnp <- subset_samples(data.filt.rarified, location == "Glacier National Park")

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
quantile(environmental_data$sodium, c(0:3/3))
                 

indval_filter <- function(in_data, target_data, cutoff = 0.05){
  ##Feature selection based on bioindicator analysis as described by 
  # Dufrene & Legendre, 1997. The target variable is split into tertiles by value
  # and the taxon / OTU table is preprocessed using Hellinger transformation
  
  dummy = target_data
  low_pos = which(dummy < 2.296667) #lower prob. must update for every env variable
  mid_pos = which(dummy >= 2.296667) #lower prob. must update for everyb env variable
  high_pos = which(dummy > 3.546667) #higher prob. must update for every env variable
  
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
write.csv(data.filt.rarified.mh@otu_table,'otus_mh.csv')
#Import phyloseq OTU table as an OTU table/dataframe
otu <-read.csv('otus_mh.csv')
#do some shuffling of the OTU table
otu_flip <- as.data.frame(t(otu)) #makes it a dataframe and puts x into y and y into x (flips it)
names(otu_flip) <- as.matrix(otu_flip[1, ]) # renames columns
otu_flip<- otu_flip[-1, ] #removes first row
otu_flip_num<-as.data.frame(lapply(otu_flip, as.numeric)) #convert from character to number
otu_flip_num$ncbi_name_18s<-row.names(otu_flip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#Repeat for metadata

#STOP! Make sure you are using the correct rarified dataset (MH, GNP, Both)
write.csv(data.filt.rarified.mh@sam_data, 'meta_18s_mh.csv')
meta <- read.csv('meta_18s_mh.csv')

## Join based on sample_id
otu_final<-left_join(meta, otu_flip_num , by = c("ncbi_name_18s" = "ncbi_name_18s")) # join based on sample IDs, assuming they're the same for both OTU table and metadata
tax_data <- otu_final[,69:5933] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
hellinger_tax_data <- as.data.frame(t(apply(tax_data, 1, function(x) sqrt(x / sum(x)))))


environmental_data <- select(otu_final, fi_best, elevation_m, ysi_temp, ysi_do, ysi_spc,
                             ysi_ph, ysi_orp, alkalinity_caco3, bicarbonate_hco3,
                             calcium, chloride, magnesium, potassium, sio2, sodium,
                             strontium, sulfate, tds)
  #the metadata column group you care about

modelfolder <- "/Users/chrishansen/Desktop/mod"
predfolder <- "/Users/chrishansen/Desktop/pred" 
evalfolder <- "/Users/chrishansen/Desktop/eval1"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "sodium", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 5,
                            feature_selection_method = "indVal")

modelfolder <- "/Users/chrishansen/Desktop/mod"
predfolder <- "/Users/chrishansen/Desktop/pred" 
evalfolder <- "/Users/chrishansen/Desktop/eval2"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "sodium", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 5,
                            feature_selection_method = "indVal")

modelfolder <- "/Users/chrishansen/Desktop/mod"
predfolder <- "/Users/chrishansen/Desktop/pred" 
evalfolder <- "/Users/chrishansen/Desktop/eval3"

##num_folds = 10 for 10-fold cross-validation; this should be the number of 
#   samples for leave-one-out cross-validation
train_covariation_framework(in_dataset = tax_data, out_dataset = environmental_data, 
                            target_feature = "sodium", model = "rf", 
                            model_target = modelfolder, prediction_target = predfolder, 
                            evaluation_target = evalfolder, num_folds = 5,
                            feature_selection_method = "indVal")


#multipatt exploration
environmental_tert <- environmental_data %>%
  mutate(recharge_level = cut(recharge_model,
                              quantile(recharge_model, c(0:3/3)),
                              labels = c("low", "medium", "high"),
                              include.lowest = TRUE))

recharge_data <- environmental_tert$recharge_level
multi_test <- multipatt(x=hellinger_tax_data, cluster=recharge_data, control = how(nperm=999))

summary(multi_test)
#correct p-values?: https://stackoverflow.com/questions/55271042/r-indicspecies-package-multipatt-function-extract-values-from-summary-multip
#i want to know which species corresponds to OTUs
multipatt_output <- multi_test$sign
mp_filter <- filter(multipatt_output, p.value <= 0.05)
#need tax data as a dataframe
write.csv(filter.data@tax_table,'tax.csv')
#Import phyloseq OTU table as an OTU table/dataframe
tax <-read.csv('tax.csv')
mp_filter$X <- rownames(mp_filter) 

join <- left_join(mp_filter, tax)
join_filter <- filter(join, s.high == 1 & s.low == 0 & s.medium == 0)
#Summarize number of bioindicators at various taxaonomic levels
join_filter %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  #top_n(n=7,wt = n) %>%
  arrange(Genus, desc(n))

#####Create rotated geom_jitter to portray R^2 variation
cov_data <- read.csv("springs_cov_rare.csv")
cov_data_regional <- filter(cov_data, cat == "Glacier National Park" | cat == "Mount Hood")

cov_data_order <- cov_data_regional %>%
  mutate(env_variable = fct_reorder(env_variable, r_squared, mean, na.rm = TRUE,
                                    .desc = TRUE))

#a failed attempt to make Fi pop
# a <- ifelse(cov_data_order$env_variable == "Fi", "bold", "plain")
# 
# ds <- data.frame(
#   a = colnames(cov_data_order),
#   b = cov_data_order$env_variable)
# 
# ds_face = rep("plain", nrow(ds))
# ds_face[106:108] <- "bold"


ggplot(cov_data_order, aes(x = env_variable, y = r_squared, color = cat)) +
  geom_jitter(alpha = 0.6, width = 0.1, size = 3) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  labs(x = "Enviornmental Variable", y = bquote(~R^2)) +
  stat_summary(fun.data = "mean_cl_boot", #Add mean and 05% confidence intervals
               show.legend = FALSE,
               position = position_nudge(x = .3, y = 0),
               fun.args=(conf.int=0.95),
               size = 1,
               pch = "o") +
  geom_hline(yintercept = summarise(filter(cov_data_order, cat == "Glacier National Park"), 
                                     mean(r_squared, na.rm = TRUE)) %>% pull(), 
            lty = 2,
            color = "#667c74",
            size = 2)+
  geom_hline(yintercept = summarise(filter(cov_data_order, cat == "Mount Hood"), 
                                    mean(r_squared, na.rm = TRUE)) %>% pull(), 
             lty = 2,
             color = "#b3692b",
             size = 2)+
  coord_flip() +
  scale_y_continuous(limits = c(0,1)) +
  theme_light() +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(axis.title = element_text(size = 14), #Modify text sizes to improve figure readability
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.y = element_blank())

