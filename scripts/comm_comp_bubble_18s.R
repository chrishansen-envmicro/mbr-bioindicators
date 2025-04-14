library(tidyverse)
library(readxl)
library(dplyr)
library(phyloseq)
library(wesanderson)
library(vegan)
library(DESeq2)
library(gridExtra)
library(GGally)
library(agricolae)
library(corncob)

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
  filter(bio_rep != "NA")
map_filt_mut <- mutate(map_filt, primary_source = 
                         case_when(fi_best > 0.6 ~ "Glacier Meltwater",
                                   fi_best <= 0.6 ~ "Precipitation"))

map_filt_mut2 <- filter(map_filt_mut, bio_rep != "blank")
map_filt_mut2 <- filter(map_filt_mut2, ncbi_name_18s != "")
map_filt_mut2 <- map_filt_mut2 %>%
  mutate(ncbi_name_16s = fct_reorder(ncbi_name_18s, primary_source))
map_filt_mut2 <- sample_data(map_filt_mut2)
rownames(map_filt_mut2) <- map_filt_mut2$ncbi_name_18s


#Merge mothur data with metadata to create phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map_filt_mut2)

#Replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Remove non-bacterial mitochondria and chloropast seq. should be removed during amp processing (mother), but this is a double-check
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

data_filter <- filter_taxa(data, function(x) sum(x)>0, TRUE) #Drop 0-count taxa
data_filt <- clean_taxa_names(data_filter) #Removes excess 0s from OTUs
#remove singletons
data_filt <- filter_taxa(data, function (x) {sum(x > 0) > 1}, prune=TRUE)

#identify which samples had fewer than 10^3 reads and remove (12)
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (6000))
min(sample_sums(data_filt_clean)) #check min at 11000. revise above to be 90% of this min (11919)

#Rarefy per alpha diversity procedure
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=0.9*min(sample_sums(data_filt_clean)), 
                                       replace=F)

#data_filt_rarified_filt <- filter_taxa(data_filt_rarified, function (x) {sum(x > 2) > 30}, prune=TRUE)

#Tranform to relative abundance
rel_abund_rarefied <- data_filt_rarified %>% transform_sample_counts(function(x) {
  x/sum(x)
})


#CPCOLS <- c("#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
          #  "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           # "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
            #"#8F5939", "#57666E", "#283C47", "#8F386F", "#660F45")

#rel_abund_bubble <- rel_abund_rarefied %>%
 # psmelt() %>%
  #filter(Abundance > 0.001) %>% #As with all filtering, you have to make and defend your choices.
  #arrange(quadrant)
ps_order <- phyloseq::tax_glom(rel_abund_rarefied, "Class")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Class"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.04, relabund = TRUE)
rel_abund_bubble_2 <- psmelt(ps_order_trim)

library(hrbrthemes)
library(viridis)
# Most basic bubble plot
rel_abund_bubble_2 %>%
  #arrange(quadrant) %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=ncbi_name_16s, y=Class, size=Abundance, fill=primary_source)) +
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(0, 24), name="Relative Abundance") +
  scale_fill_manual(values=wes_palette("IsleofDogs1"), name = "Sourcing") +
  facet_wrap(~location, nrow = 1) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2), axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18), legend.text =  element_text(size = 15), legend.title =  element_text(size = 15))+
  ylab("Order") +
  xlab("Sample") +
  theme(axis.title = element_text(size = 14))








#We first agglomerate data at Genus level and then create a bar plot object
ps_order = tax_glom(rel_abund_rarefied, "Order")
ps_bar_order <- plot_bar(ps_order) 

ps_bar_order$layers <- ps_bar_order$layers[-1]

# sets maximum relative abundance for figure legend dimensions
maxC<-max(otu_table(ps_order)) 

# replace with bubble plot instead of bars
ps_bubble_order <- ps_bar_order + 
  geom_point(aes_string(x="alt_sample_name", y="Order", size="Abundance", color = "quadrant"), alpha = 0.7) +
  scale_size_continuous(limits = c(0.001,maxC)) +
  xlab("Sample") +
  ylab("Order") +
  ggtitle("Relative abundances at Genus level") + 
  labs(caption = "Abundances below 0.001 were considered absent") +
  theme_light()

# Plot
ps_bubble_order

newtab = data.table(ps_bubble_order$data)
setorder(newtab, quadrant)
ps_bubble_order$data <- newtab
ps_bubble_order

library("data.table")
newtab = data.table(pTime$data)
setorder(newtab, TIMEPOINT_NUMBER)
pTime$data <- newtab
print(pTime)

ggplot(ps_melt_ord_ord, aes(x = factor(Sample), y = Abundance*100, fill = Order)) +
  geom_bar(stat = "identity") +
  facet_wrap(~primary_source, scales = "free_x") +
  theme_light() +
  #scale_fill_manual(values=CPCOLS) +
  theme(axis.title.x= element_blank(), axis.text.x=element_text(angle=90, hjust=1))+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Percent Relative Abundance") +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        legend.position = "bottom", legend.direction = "horizontal")

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


#Filter OTUs that were at at least a ___ relative abundance within samples
data_filt_rarefied_trim <- phyloseq_filter_sample_wise_abund_trim(data_filt_rarified, minabund = 0.005, relabund = TRUE)
melt_1 <- psmelt(data_filt_rarefied_trim)

#Create boolean table (TRUE = >1% relative abundance)
rel_abund_rarefied_bool = filter_taxa(rel_abund_rarefied, function(x) (sum(x)/45) >= 0.0000001, FALSE)

#Drop taxa with <1% relative abundance
rel_abund_rarefied_filt = prune_taxa(rel_abund_rarefied_bool, rel_abund_rarefied)

###### Agglomerate to phylum-level and rename ######
ps_phylum <- phyloseq::tax_glom(rel_abund_rarefied, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
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

#Count of phyla that compose 90% of the relative abundance
a[1:15,] %>%
  summarise(sum = sum(mean))

#Check to see if Phyla differ by recharge source
ps_melt_phylum %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  #scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr"))

#Check to see if Phyla differ by location
ps_melt_1 %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH"))

##### Agglomerate to Class-level and rename ######
ps_class <- phyloseq::tax_glom(rel_abund_rarefied, "Class")
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

#Check to see if Phyla differ by recharge source
ps_melt_class %>%
  ggplot(aes(x = primary_source, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier Meltwater" = "fi",
                            "Precipitation" = "fs/fr"))

#Check to see if Phyla differ by location
ps_melt_class %>%
  ggplot(aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Abundance\n", color = "Phyla") +
  facet_wrap(~ OTU, scales = "free") +
  scale_x_discrete(labels=c("Glacier National Park" = "GNP",
                            "Mount Hood" = "MH"))

##### Agglomerate to Order-level and rename ######
ps_order <- phyloseq::tax_glom(rel_abund_rarefied, "Order")
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
