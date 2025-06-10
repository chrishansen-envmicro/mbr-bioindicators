#load libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(wesanderson)
library(ggpubr)
library(data.table)
library(cowplot)
library(rstatix)
library(ggfortify)
library(ggtext)
library(broom)


# reference(s)
# Miller, J. B., Frisbee, M. D., Hamilton, T. L., & Murugapiran, S. K. (2021). 
# Recharge from glacial meltwater is critical for alpine springs and their microbiomes. 
# Environmental Research Letters, 16(6), 64012-. https://doi.org/10.1088/1748-9326/abf06b 

#code sourcing (partial): https://rpubs.com/lconteville/713954

# import mothur output produced from Miller et al. (2021)
sharedfile <- "16S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "16S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated.csv"

# returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
# import metadata
map <- read.csv(mapfile)

# rename location values to improve figure clarity
map$location <- recode(map$location, MH = "Mount Hood") %>% 
  recode(GNP = "Glacier National Park") 

#filter to extract variables of interest and avoid duplicate env values from replicates
map_filt <- filter(map, sample_type == "Spring" | sample_type == "Stream") %>%
  filter(ncbi_name_16s != "")
map_filt_mut <- mutate(map_filt, primary_source = 
                         case_when(fi_best > 0.6 ~ "Glacier Meltwater",
                                   fi_best <= 0.6 ~ "Precipitation"))

map_filt_mut <- sample_data(map_filt_mut)
rownames(map_filt_mut) <- map_filt_mut$ncbi_name_16s


#Merge mothur data with metadata to create phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map_filt_mut)

#Replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#Remove non-bacterial mitochondria and chloropast seq. should be removed during amp processing (mother), but this is a double-check
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

data_filter <- filter_taxa(data, function(x) sum(x)>0, TRUE) #Drop 0-count taxa
data_filt <- clean_taxa_names(data_filter) #Removes excess 0s from OTUs
#remove singletons
data_filt <- filter_taxa(data, function (x) {sum(x > 0) > 1}, prune=TRUE)

#identify which samples had fewer than 10^3 reads and remove (12)
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (7000))
min(sample_sums(data_filt_clean)) #check min at 11000. revise above to be 90% of this min (11919)

#Rarefy per alpha diversity procedure
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=0.9*min(sample_sums(data_filt_clean)), 
                                       replace=F)

#Tranform to relative abundance - rare
rel_abund_rarefied <- data_filt_rarified %>% transform_sample_counts(function(x) {
  x/sum(x)
})

#Tranform to relative abundance - no rare
rel_abund_clean <- data_filt_clean %>% transform_sample_counts(function(x) {
  x/sum(x)
})

#Use psmelt to support formatting demands
#rel_abund_data_melt <- psmelt(rel_abund_data)
#write.csv(rel_abund_data_melt,'melt.csv')

#Summary statistics for community data
#Number of ...
#unique(melt_1$Phylum) 
#unique(rel_abund_data_melt$Class)
#unique(rel_abund_data_melt$Order)
#unique(rel_abund_data_melt$Family)
#unique(rel_abund_data_melt$Genus)

#Quickly see differences in relative abundance of Phyla
#https://github.com/joey711/phyloseq/issues/1501
#https://github.com/joey711/phyloseq/issues/1129
#https://github.com/joey711/phyloseq/issues/1089#issuecomment-471334036
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html

###### Agglomerate to phylum-level and rename ######
ps_phylum <- phyloseq::tax_glom(rel_abund_rarefied, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
# from "comm_fxns.R" script
ps_phylum_trim <- phyloseq_filter_sample_wise_abund_trim(ps_phylum, minabund = 0.02, relabund = TRUE)

#Melt and plot
ps_melt_phylum <- phyloseq::psmelt(ps_phylum_trim)

#Summary statistics for community composition
ps_melt_phylum %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

CPCOLS <- c("#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
            "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
            "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
            "#8F5939", "#8F386F", "#660F45", "#57666E", "#283C47") 
#Check to see if Phyla differ by recharge source
ps_melt_phylum %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr")) +
  theme_light()

#Check to see if Phyla differ by location
ps_melt_phylum %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH"))+
  theme_light()

##### Agglomerate to Class-level and rename ######
ps_class <- phyloseq::tax_glom(rel_abund_clean, "Class")
ps_class_edit <- subset_taxa(ps_class, Class != "uncultured")
phyloseq::taxa_names(ps_class_edit) <- phyloseq::tax_table(ps_class_edit)[,"Class"]
phyloseq::otu_table(ps_class_edit)[1:5, 1:5]
ps_class_trim <- phyloseq_filter_sample_wise_abund_trim(ps_class_edit, minabund = 0.04, relabund = TRUE)


#Melt and plot
ps_melt_class <- phyloseq::psmelt(ps_class_trim)

#Summary statistics for community composition
ps_melt_class %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

#Check to see if Classes differ by recharge source
ps_melt_class %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr")) +
  theme_light()

#Check to see if Classes differ by location
ps_melt_class %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH")) +
  theme_light()


##### Agglomerate to Order-level and rename ######
ps_order <- phyloseq::tax_glom(rel_abund_clean, "Order")
ps_order_edit <- subset_taxa(ps_order, Order != "uncultured" & Class != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Order"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.06, relabund = TRUE)


#Melt and plot
ps_melt_order <- phyloseq::psmelt(ps_order_trim)

#Summary statistics for community composition
ps_melt_order %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

#Check to see if Orders differ by recharge source
ps_melt_order %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr")) +
  theme_light()

#Check to see if Phyla differ by location
ps_melt_order %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH")) +
  theme_light()

##### Agglomerate to Family-level and rename ######
ps_family<- phyloseq::tax_glom(rel_abund_rarefied, "Family")
ps_family_edit <- subset_taxa(ps_family, Order != "uncultured" & 
                                Class != "uncultured" &
                                Family != "uncultured" &
                                Family != "Unknown_Family")
phyloseq::taxa_names(ps_family_edit) <- phyloseq::tax_table(ps_family_edit)[,"Family"]
phyloseq::otu_table(ps_family_edit)[1:5, 1:5]
ps_family_trim <- phyloseq_filter_sample_wise_abund_trim(ps_family_edit, 
                                                         minabund = 0.06, relabund = TRUE)


#Melt and plot
ps_melt_order <- phyloseq::psmelt(ps_family_trim)

#Summary statistics for community composition
ps_melt_order %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

#Check to see if Phyla differ by recharge source
ps_melt_order %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr"))

#Check to see if Phyla differ by location
ps_melt_order %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH"))

##### Agglomerate to Genus-level and rename ######
ps_family<- phyloseq::tax_glom(rel_abund_rarefied, "Genus")
ps_family_edit <- subset_taxa(ps_family, Order != "uncultured" & 
                                Class != "uncultured" &
                                Family != "uncultured" &
                                Genus != "uncultured")
phyloseq::taxa_names(ps_family_edit) <- phyloseq::tax_table(ps_family_edit)[,"Genus"]
phyloseq::otu_table(ps_family_edit)[1:10, 1:10]
ps_family_trim <- phyloseq_filter_sample_wise_abund_trim(ps_family_edit, 
                                                         minabund = 0.06, relabund = TRUE)


#Melt and plot
ps_melt_order <- phyloseq::psmelt(ps_family_trim)

#Summary statistics for community composition
ps_melt_order %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

#Check to see if Phyla differ by recharge source
ps_melt_order %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr"))

#Check to see if Phyla differ by location
ps_melt_order %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH"))

####Summary table
library(kableExtra)
abund <- read_csv("abund_stats.csv")

kbl(abund) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F, font_size = 16) %>%
  column_spec(1:2, bold = T) %>%
  row_spec(1, bold = T, color = "white", background = "#4092A8") %>%
  row_spec(2, bold = T, color = "white", background = "#549EB1") %>%
  row_spec(3, bold = T, color = "white", background = "#67A9B5") %>%
  row_spec(6, bold = T, color = "white", background = "#4092A8") %>%
  row_spec(7, bold = T, color = "white", background = "#549EB1") %>%
  row_spec(8, bold = T, color = "white", background = "#67A9B5") %>%
  row_spec(9, bold = T, color = "white", background = "#4092A8") %>%
  row_spec(10, bold = T, color = "white", background = "#549EB1") %>%
  row_spec(11, bold = T, color = "white", background = "#67A9B5") %>%
  row_spec(12, bold = T, color = "white", background = "#4092A8") %>%
  row_spec(13, bold = T, color = "white", background = "#549EB1") %>%
  row_spec(14, bold = T, color = "white", background = "#67A9B5") %>%
  row_spec(16, bold = T, color = "white", background = "#4092A8") %>%
  row_spec(17, bold = T, color = "white", background = "#549EB1") %>%
  row_spec(18, bold = T, color = "white", background = "#67A9B5")

#Exploratory venndiagrams
summary_plot_phyloseq(data_filt_rarified)
otu_circle(data_filt_rarified, 'location', taxa = "Genus")

#Attempt2
#https://microbiome.github.io/tutorials/core_venn.html
#convert to relative abundances
library(microbiome)
library(microbiomeutilities)
library(eulerr)
pseq.rel <- microbiome::transform(data_filt_rarified, "compositional")

#make a list of variables i am comparing
locs <- unique(as.character(meta(pseq.rel)$primary_source))
print(locs)

#Write a for loop to go through each of the disease_states one by one and 
#combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in locs){ # for each variable n in locs
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, primary_source == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, # 0.001 in atleast 90% samples 
                         prevalence = 0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
#print venn
mycols <- c("Precipitation"="#d8a499", "Glacial Ice"="#c5cdf7") 
plot(venn(list_core),
     fill = mycols)

#Attempt2 - LOCATION
#https://microbiome.github.io/tutorials/core_venn.html
#convert to relative abundances
pseq.rel <- microbiome::transform(data_filt_rarified, "compositional")

#make a list of variables i am comparing
locs <- unique(as.character(meta(pseq.rel)$location))
print(locs)

#Write a for loop to go through each of the disease_states one by one and 
#combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in locs){ # for each variable n in locs
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, location == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, # 0.001 in atleast 90% samples 
                         prevalence = 0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
#print venn
mycols <- c("Glacier National Park"="#b3692b", "Mount Hood"="#667c74") 
plot(venn(list_core),
     fill = mycols)
