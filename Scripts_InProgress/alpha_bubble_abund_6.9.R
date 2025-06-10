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
mapfile <- "Master_Data_Sheet_11112021_consolidated_4.12.csv"

# returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
# import metadata
map <- read.csv(mapfile)

# rename location values to improve figure clarity
map$location <- recode(map$location, MH = "Mount Hood") %>% 
  recode(GNP = "Glacier National Park") 

##### for relative abundance plots ######
#filter to extract variables of interest and avoid duplicate env values from replicates
map_filt <- filter(map, ncbi_name_16s != "")
map_filt_mut <- mutate(map_filt, primary_source = 
                         case_when(fi_best > 0.6 ~ "Glacier Meltwater",
                                   fi_best <= 0.6 ~ "Precipitation"))

map_filt_mut <- sample_data(map_filt_mut)
rownames(map_filt_mut) <- map_filt_mut$ncbi_name_16s

#filter to extract variables of interest and avoid duplicate env values from replicates
map_filt <- filter(map, ncbi_name_16s != "")
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

library(corncob)
data_filt <- clean_taxa_names(data_filter) #Removes excess 0s from OTUs

# remove singletons
data_filt <- filter_taxa(data, function (x) {sum(x > 0) > 1}, prune=TRUE)

# identify which samples had fewer than 8000 reads and remove - See methods.
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (7000))
min(sample_sums(data_filt_clean)) 

# rarefy per alpha diversity procedure
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=min(sample_sums(data_filt_clean)), 
                                       replace=F)

data_filt_rarified_mh = subset_samples(data_filt_rarified, location == "Mount Hood")
data_filt_rarified_gnp = subset_samples(data_filt_rarified, location == "Glacier National Park")

# how many samples remain from each mountain range? --> threshold created to align with alpha-diversity sample size
sample_data(data_filt_rarified) %>%
  group_by(location) %>%
  summarise(n = n())

# REPEAT for only springs and streams
stream.spring <- c("Stream", "Spring")
data_spst <- subset_samples(data_filt_clean, sample_type %in% stream.spring)

data_filt_rarified_spst = rarefy_even_depth(data_spst, rngseed=1, 
                                            sample.size=min(sample_sums(data_filt_clean)), 
                                            replace=F)

data_filt_rarified_mh_spst = subset_samples(data_filt_rarified_spst, location == "Mount Hood")
data_filt_rarified_gnp_spst = subset_samples(data_filt_rarified_spst, location == "Glacier National Park")

# how many samples remain from each mountain range? --> threshold created to align with alpha-diversity sample size
sample_data(data_filt_rarified_spst) %>%
  group_by(location) %>%
  summarise(n = n())


## tranform to relative abundance

# all samples
rel_abund_rarefied <- data_filt_rarified_spst %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# MH only
rel_abund_rarefied_mh <- data_filt_rarified_mh_spst %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# GNP only
rel_abund_rarefied_gnp <- data_filt_rarified_mh_spst %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# total - no rare
rel_abund_norarefied <- data_spst %>% transform_sample_counts(function(x) {
  x/sum(x)
})


# MH only - no rare
data_filt_clean_mh <- subset_samples(data_spst, location == "Mount Hood")
rel_abund_norarefied_mh <- data_filt_clean_mh %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# GNP only - no rare 
data_filt_clean_gnp <- subset_samples(data_spst, location == "Glacier National Park")
rel_abund_norarefied_gnp <- data_filt_clean_gnp %>% transform_sample_counts(function(x) {
  x/sum(x)
})


# PS Melt to consolidate abundance at specified taxanomic level (e.g., Order)

# both systems
ps_order <- phyloseq::tax_glom(rel_abund_norarefied, "Phylum")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Phylum"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.004, relabund = TRUE)
rel_abund_bubble_2 <- psmelt(ps_order_edit)

# MH

# Extract Family names
family_names <- as.vector(phyloseq::tax_table(ps_order_edit)[, "Family"])

# Check for duplicates
duplicated_families <- family_names[duplicated(family_names)]

ps_order <- phyloseq::tax_glom(rel_abund_norarefied_mh, "Phylum")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Phylum"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
rel_abund_bubble_2_mh <- psmelt(ps_order_edit)

# summary statistics for community composition
a <- rel_abund_bubble_2_mh %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance)*100,
            min = min(Abundance)*100,
            max = max(Abundance)*100) %>%
  arrange(desc(mean))

# GNP
ps_order <- phyloseq::tax_glom(rel_abund_norarefied_mh, "Genus")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured"
                             & Order != "uncultured" & Family != "uncultured"
                             & Family != "Unknown_Family" & Genus != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Genus"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
rel_abund_bubble_2_mh <- psmelt(ps_order_edit)

# summary statistics for community composition
b <- rel_abund_bubble_2_mh %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance)*100,
            min = min(Abundance)*100,
            max = max(Abundance)*100) %>%
  arrange(desc(mean))

# determine how many taxa comprise 90% of the total relative abundance
sum(b$mean[1:10])


rel_abund_norarefied

# for MH
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_mh, "Phylum")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Order"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
#ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
#                                                        minabund = 0.045, relabund = TRUE)

rel_abund_bubble_mh <- psmelt(ps_order_trim)

#for GNP
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_gnp, "Genus")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured"
                             & Class != "NA" & Order != "uncultured" & Family != "uncultured"
                             & Family != "Unknown_Family" & Genus != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Genus"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.004, relabund = TRUE)

rel_abund_bubble_gnp <- psmelt(ps_order_trim)


a <- as.data.frame(tax_table(ps_order_edit))
duplicated(a$Family)
a$Family[76]

# summary statistics for community composition
rel_abund_bubble_gnp %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))

#Produce bubble plot
CPCOLS <- c("#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
            "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
            "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
            "#8F5939", "#8F386F", "#660F45", "#57666E", "#283C47")

# generate relative abundance bar graph for Figure S10.
ggplot(rel_abund_bubble_2, aes(x = factor(alt_sample_name), y = Abundance*100, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~location, scales = "free_x") +
  theme_light() +
  scale_fill_manual(values=CPCOLS) +
  theme(axis.title.x= element_blank(), axis.text.x=element_text(angle=90, hjust=1))+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Percent Relative Abundance") +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        legend.position = "bottom", legend.direction = "horizontal",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.text = element_text(size = 13))



library(hrbrthemes)
library(viridis)
library(see)
# Most basic bubble plot
write.csv(rel_abund_bubble_2, file = "bubble.csv")
bubble <- read.csv("bubble.csv")

rel_abund_bubble_2 %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=alt_sample_name, y=Class, size=Abundance, shape = primary_source)) +
  geom_point(alpha=0.7, aes(shape = primary_source, color = sample_type)) +
  scale_size(range = c(0.1, 15), name="Relative Abundance") +
  scale_fill_okabeito(palette = "full") +
  theme_light() +
  ylab("Order") +
  xlab("Sample") +
  theme_ipsum() +
  labs(color = "Sample Type", shape = "Primary Source") +
  #facet_wrap(~location) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.x= element_text(size = 18),
        axis.title.y= element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 12),
        axis.text.x=element_text(angle=90, hjust=1),
        legend.key.size = unit(1,"line"),
        legend.position = "right")

desired_sample_type_order <- c("Glacial Melt", "Snow", "Snow Algae", "Spring", "Stream")

library(forcats)
library(dplyr)

#formatting (mh)
rel_abund_bubble_2_mh <- rel_abund_bubble_2_mh %>%
  mutate(
    sample_type = factor(sample_type, levels = desired_sample_type_order)
  ) %>%
  arrange(sample_type, alt_sample_name) %>%
  mutate(
    alt_sample_name = factor(alt_sample_name, levels = unique(alt_sample_name))
  ) %>%
  filter(Abundance > 0)

#formatting (gnp)
rel_abund_bubble_2_gnp <- rel_abund_bubble_2_gnp %>%
  mutate(
    sample_type = factor(sample_type, levels = desired_sample_type_order)
  ) %>%
  arrange(sample_type, alt_sample_name) %>%
  mutate(
    alt_sample_name = factor(alt_sample_name, levels = unique(alt_sample_name))
  ) %>%
  filter(Abundance > 0)

# Let's embolden discussed Classes
# install.packages("ggtext") if not already installed
library(ggtext)

# Let's say you want to bold "Gammaproteobacteria" and "Bacteroidia"
rel_abund_bubble_2_mh$Class <- factor(rel_abund_bubble_2_mh$Class)

# Get unique levels
class_levels <- levels(rel_abund_bubble_2_mh$Class)

# Define bold labels
custom_labels <- setNames(class_levels, class_levels)
bold_these <- c("Gammaproteobacteria", "Bacteroidia", "Anaerolineae", "Nitrososphaeria",
                "Nanoarchaeia", "Alphaproteobacteria")
custom_labels[bold_these] <- paste0("<b>", bold_these, "</b>")

#Primary Plot (Fig 7)

rel_abund_bubble_2_mh %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=alt_sample_name, y=Class, size=Abundance, shape = primary_source)) +
  geom_point(alpha=0.7, aes(shape = primary_source, color = sample_type), stroke = 1.5) +
  scale_size(range = c(0.1, 15), name="Relative Abundance") +
  scale_color_manual(values = c("#85b6b2","#7b6e97", "#7294d4", "#a2a475", "#c44e52")) +
  scale_shape_manual(values = c(19,1))+
  theme_light() +
  ylab("") +
  xlab("") +
  # theme_ipsum() +
  labs(color = "Sample Type", shape = "Primary Source") +
  #facet_wrap(~location) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.x= element_text(size = 22),
        axis.title.y= element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 13),
        axis.text.x=element_text(angle=90, hjust=1, size = 14),
        axis.text.y = ggtext::element_markdown(size = 14),
        legend.key.size = unit(1,"line"),
        legend.position = "right") +
  scale_y_discrete(labels = custom_labels)


rel_abund_bubble_2_gnp %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=alt_sample_name, y=Class, size=Abundance, shape = primary_source)) +
  geom_point(alpha=0.7, aes(shape = primary_source, color = sample_type), stroke = 1.5) +
  scale_size(range = c(0.1, 15), name="Relative Abundance") +
  scale_color_manual(values = c("#85b6b2", "#a2a475", "#c44e52")) +
  scale_shape_manual(values = c(19,1))+
  theme_light() +
  ylab("") +
  xlab("") +
  theme_ipsum() +
  labs(color = "Sample Type", shape = "Primary Source") +
  #facet_wrap(~location) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.x= element_text(size = 22),
        axis.title.y= element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 13),
        axis.text.x=element_text(angle=90, hjust=1, size = 14),
        axis.text.y=element_text(size = 14),
        legend.key.size = unit(1,"line"),
        legend.position = "right")

rel_abund_bubble_2 %>%
  ggplot(aes(x=alt_sample_name, y=Phylum, size=Abundance, shape = primary_source)) +
  geom_point(alpha=0.7, aes(shape = primary_source, color = sample_type)) +
  scale_size(range = c(0.1, 15), name="Relative Abundance") +
  scale_color_manual(values = c("#E69F00", "#F5C710", "#0072B2" )) +
  theme_light() +
  ylab("Order") +
  xlab("Sample") +
  theme_ipsum() +
  labs(color = "Sample Type", shape = "Primary Source") +
  #facet_wrap(~location) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title.x= element_text(size = 18),
        axis.title.y= element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 12),
        axis.text.x=element_text(angle=90, hjust=1),
        legend.key.size = unit(1,"line"),
        legend.position = "right")

##### LEFse Analysis #####
library(microbiomeMarker)

# mh sourcing differences
lef_out<-run_lefse(data_filt_rarified_mh, group = "primary_source", norm = "none", 
                   kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out
plot_ef_bar(lef_out)

# gnp sourcing differences
lef_out<-run_lefse(data_filt_rarified_gnp, group = "primary_source", norm = "none", 
                   kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out
plot_ef_bar(lef_out)


## universal sourcing differences
lef_out_both <-run_lefse(data_filt_rarified, group = "primary_source", norm = "none", 
                         kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out_both
plot_ef_bar(lef_out_both)

## geographic differences
lef_out_both <-run_lefse(data_filt_rarified, group = "location", norm = "none", 
                         kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out_both
plot_ef_bar(lef_out_both)

# mh quadrant differences

#remove NW quadrant as there is only 1 sample
data_filt_rarified_mh_nonw <- subset_samples(data_filt_rarified_mh, quadrant != "NW")

lef_out_reg <-run_lefse(data_filt_rarified_mh_nonw, group = "quadrant_ew", norm = "none", 
                        kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out_reg
plot_ef_bar(lef_out_reg)

# mh quadrant differences
lef_out_quad <-run_lefse(data_filt_rarified_mh_quad, group = "quadrant", norm = "none", 
                         kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out_quad
plot_ef_bar(lef_out_quad)

sample_data(data_filt_rarified_mh) %>%
  group_by(quadrant) %>%
  summarise(n = n())

sample_data(data_filt_rarified_mh) %>%
  group_by(region) %>%
  summarise(n = n())

data_filt_rarified_mh_reg <- subset_samples(data_filt_rarified_mh, region != "Eliot Region")
data_filt_rarified_mh_quad <- subset_samples(data_filt_rarified_mh, quadrant != "NW")


# Let's identify biomarkers between melt and precip end-members for use in our relative abundance plot
# dat_type<-subset_samples(data_filt_rarified, samp_type %in% c("Glacial Melt","Spring_Stream"))
# meta_table(dat_type)

lef_out_em<-run_lefse(data_filt_rarified, group = "location", norm = "none",
                      kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out_em
plot_ef_bar(lef_out_em)

#https://microbiome.github.io/tutorials/core_venn.html
#convert to relative abundances
library(microbiome)
library(gplots)
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
                         detection = 0.0, # 0.001 in atleast 90% samples 
                         prevalence = 0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
#print venn
# Make sure the names in 'mycols' match the set names in 'list_core'
mycols <- c("Mount Hood" = "#b3692b", "Glacier Meltwater" = "#667c74")

# Generate Venn diagram
plot(eulerr::venn(list_core),
     fills = list(fill = mycols, alpha = 0.7),
     edges = list(col = "gray20", lwd = 1),
     labels = list(col = "black", font = 2, cex = 1.3),
     quantities = list(col = "black", font = 1),
     adjust_labels = TRUE)

######## Repeat by location

pseq.rel <- microbiome::transform(data_filt_clean, "compositional")

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
mycols <- c("Precipitation"="#ED96FF", "Glacier Meltwater"="#00DEE2") 
plot(venn(list_core),
     fills = mycols)

#attempt miscmetabar: https://github.com/adrientaudiere/Miscmetabar
#summary of what is in my phyloseq object
library(MiscMetabar)
summary_plot_phyloseq(filter.data.spst)

otu_circle(filter.data.spst, 'location', taxa = "Class")




