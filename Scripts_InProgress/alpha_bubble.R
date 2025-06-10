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

library(corncob)
data_filt <- clean_taxa_names(data_filter) #Removes excess 0s from OTUs

# remove singletons
data_filt <- filter_taxa(data, function (x) {sum(x > 0) > 1}, prune=TRUE)

# identify which samples had fewer than 8000 reads and remove - See methods.
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (8000))
min(sample_sums(data_filt_clean)) 

# rarefy per alpha diversity procedure
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=0.9*min(sample_sums(data_filt_clean)), 
                                       replace=F)

data_filt_rarified_mh = subset_samples(data_filt_rarified, location == "Mount Hood")
data_filt_rarified_gnp = subset_samples(data_filt_rarified, location == "Glacier National Park")

# how many samples remain from each mountain range? --> threshold created to align with alpha-diversity sample size
sample_data(data_filt_rarified) %>%
  group_by(location) %>%
  summarise(n = n())

## tranform to relative abundance

# all samples
rel_abund_rarefied <- data_filt_rarified %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# MH only
rel_abund_rarefied_mh <- data_filt_rarified_mh %>% transform_sample_counts(function(x) {
  x/sum(x)
})

# GNP only
rel_abund_rarefied_gnp <- data_filt_rarified_gnp %>% transform_sample_counts(function(x) {
  x/sum(x)
})




# PS Melt to consolidate abundance at specified taxanomic level (e.g., Order)

# both systems
ps_order <- phyloseq::tax_glom(rel_abund_rarefied, "Order")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured"
                             & Order != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Order"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.025, relabund = TRUE)
rel_abund_bubble_2 <- psmelt(ps_order_trim)

# MH
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_mh, "Class")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Class"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.02, relabund = TRUE)
rel_abund_bubble_2_mh <- psmelt(ps_order_trim)

# GNP
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_gnp, "Class")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured" & Class != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Class"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.02, relabund = TRUE)
rel_abund_bubble_2_gnp <- psmelt(ps_order_trim)

# summarize mean abundance pof 63x phyla found in our dataset
phylum_abundance <- rel_abund_bubble_2 %>%
  group_by(Phylum) %>%
  summarize(mean_abundance = mean(Abundance) * 100) %>%
  arrange(desc(mean_abundance))

# determine how many taxa comprise 90% of the total relative abundance
sum(phylum_abundance$mean_abundance[1:10])

# summary statistics for community composition
rel_abund_bubble_2 %>%
  group_by(OTU) %>%
  summarise(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance)) %>%
  arrange(desc(mean))



# for MH
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_mh, "Phylum")
ps_order_edit <- subset_taxa(ps_order, Phylum != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Order"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
#ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
#                                                        minabund = 0.045, relabund = TRUE)

rel_abund_bubble_mh <- psmelt(ps_order_trim)

#for GNP
ps_order <- phyloseq::tax_glom(rel_abund_rarefied_gnp, "Order")
ps_order_edit <- subset_taxa(ps_order, Order != "uncultured" & Class != "uncultured")
phyloseq::taxa_names(ps_order_edit) <- phyloseq::tax_table(ps_order_edit)[,"Order"]
phyloseq::otu_table(ps_order_edit)[1:5, 1:5]
ps_order_trim <- phyloseq_filter_sample_wise_abund_trim(ps_order_edit, 
                                                        minabund = 0.1, relabund = TRUE)

rel_abund_bubble_gnp <- psmelt(ps_order_trim)

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

rel_abund_bubble_2_gnp <- rel_abund_bubble_2_gnp %>%
  mutate(
    sample_type = factor(sample_type, levels = desired_sample_type_order)
  ) %>%
  arrange(sample_type, alt_sample_name) %>%
  mutate(
    alt_sample_name = factor(alt_sample_name, levels = unique(alt_sample_name))
  ) %>%
  filter(Abundance > 0)

rel_abund_bubble_2_gnp %>%
  #mutate(country = factor(country, country)) %>%
  ggplot(aes(x=alt_sample_name, y=Class, size=Abundance, shape = primary_source)) +
  geom_point(alpha=0.7, aes(shape = primary_source, color = sample_type), stroke = 1.5) +
  scale_size(range = c(0.1, 15), name="Relative Abundance") +
  scale_color_manual(values = c("#85b6b2","#7b6e97", "#7294d4", "#a2a475", "#c44e52")) +
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

# mh regional differences
lef_out_reg <-run_lefse(data_filt_rarified_mh_reg, group = "region", norm = "none", 
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


