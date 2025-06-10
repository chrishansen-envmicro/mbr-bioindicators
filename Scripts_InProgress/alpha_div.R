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

map.recharge <- filter(map, ncbi_name_16s != "") # remove samples without 16S data
map.recharge <- sample_data(map.recharge) # format sample data
rownames(map.recharge) <- map.recharge$ncbi_name_16s # align rownames across metadata .tax and .shared files

# merge morthur files (.taxonomy, .shared) with metadata
moth_merge <- merge_phyloseq(mothur_data, map.recharge)

# replace rank with taxonomic categories
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# remove non-bacterial mitochondria and chloroplast seqs
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

sum(otu_table(data)) # n_reads
nsamples(data) # n_samples
ntaxa(data) # n_OTUs

# filter data to springs and streams from MH and GNP
stream.spring <- c("Stream", "Spring")
data_spst <- subset_samples(data, sample_type %in% stream.spring)

# remove singletons
data_filt <- filter_taxa(data_spst, function (x) {sum(x > 0) > 1}, prune=TRUE)

# develop rarefaction curves and visualize distribution of read lengths
cols <-  c("#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F", "#660F45", "#707070",
           "#333333", "black", "#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F", "#660F45", "#707070",
           "#333333", "black", "#EDB4B4", "#F56C6C", "#F0C287", "#E87E37", "#F0E68C",
           "#E0DB3E", "#9DC19F", "#538257", "#B7E8E8", "#729E9D", "#CFBCD6",
           "#80588A", "#F2B8DF", "#AB1B7B", "#C79797", "#995050", "#AB9A74",
           "#8F5939", "#57666E", "#283C47", "#8F386F")

line <- c(rep(1,25), rep(4,25), rep(6,21)) # vary line type to avoid excessive colors

otu_tab <- otu_table(data_filt) # extract otu table
class(otu_tab) <- "matrix" # as.matrix() will do nothing. you get a warning here, but this is what we need to have
otu_tab <- t(otu_tab) # transpose observations to rows
rare <- rarecurve(otu_tab, step=100, ylab="OTU",  label=T, col = cols, lty = line)
abline(v=7000, lw = 1) # display selected cutoff across samples that balance 

# how many OTUs in each sample?
df <- otu_tab %>% 
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")
otu_summary <- df %>%
  group_by(SampleID) %>%
  summarize(non_zero_otus = sum(across(where(is.numeric), ~ .x != 0)))

# how many reads in each sample?
row_sums <- rowSums(otu_tab)
sorted_row_sums <- sort(row_sums)
print(sorted_row_sums)

# identify which samples had fewer than 7000 reads and remove
data_filt_clean <- subset_samples(data_filt, sample_sums(data_filt) >= (7000))
min(sample_sums(data_filt_clean)) # check min at 7000
sort(sample_sums(data_filt))

# how many samples remain from each mountain range?
sample_data(data_filt_clean) %>%
  group_by(location) %>%
  summarise(n = n())

# rarefy the samples without replacement. Rarefaction is used to simulate even number 
# of reads per sample. In this example, the rarefaction depth chosen is the 90% of the 
# minimum sample depth in the dataset (in this case 459 reads per sample).
# rarefy based on lowers read count about 1k. (1015)
# https://micca.readthedocs.io/en/latest/phyloseq.html
data_filt_rarified = rarefy_even_depth(data_filt_clean, rngseed=1, 
                                       sample.size=min(sample_sums(data_filt_clean)), 
                                       replace=F)

sum(otu_table(data_filt_rarified)) # n_reads
nsamples(data_filt_rarified) # n_samples
ntaxa(data_filt_rarified) # n_OTUs

# provides alpha diversity estimates for each sample
alphadiv <- estimate_richness(data_filt_rarified)

## test assumptions for statistical testing (i.e., parametric vs. non-parametric data?)
alphadiv_form <- rownames_to_column(alphadiv, var = "ncbi_name_16s") # extract rownames as a column for left_join
alphadiv_form <- alphadiv_form %>% # "VLOOKUP" to add location
  left_join(map %>% select(ncbi_name_16s, location), by = "ncbi_name_16s")

# create linear models for hypotheses tested to assess normality
alpha_lm_rich <- lm(Observed ~ location, data = alphadiv_form)
alpha_lm_shan <- lm(Shannon ~ location, data = alphadiv_form)
autoplot(alpha_lm_rich) # Near-normal, equal variance of residuals --> t-test
autoplot(alpha_lm_shan) # deviates from normal, heteroscedasticity --> wilcox rank-sum test

# test the hypotheses that mean alpha diversity metrics differ between sites
t.test(Observed ~ location, data = alphadiv_form, var.equal = TRUE)
wilcox.test(Shannon ~ location, data = alphadiv_form)

# visualizing mean differences for metrics and sd differences
alphadiv_form %>%
  group_by(location) %>%
  summarize(sd_shann = sd(Shannon),
            mean_shann = mean(Shannon),
            sd_obs = sd(Observed),
            mean_obs = mean(Observed))


# summarize standard deviation
alphadiv_form %>%
  group_by(location) %>%
  summarise(observed_sd = sd(Observed),
            shannon_sd = sd(Shannon))

# define comparison for statistical analysis
my_comparisons <- list( c("Glacier National Park", "Mount Hood")) 
# manually add Wilcox values for the hypothesis-testing aesthetics
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# plot alpha diversity using Observed (richness) and Shannon (Shannon diversity) metrics
p <- plot_richness(data_filt_rarified, x="location", 
                   measures=c("Observed", "Shannon"),
                   color = "location") 

p$layers <- p$layers[-1]

# test comparisons if removing replicates
data_filt_rarified_noreps <- subset_samples(data_filt_rarified,
                                            bio_rep != 2 & bio_rep != 3)
p <- plot_richness(data_filt_rarified, x="location", 
                        measures=c("Observed", "Shannon"),
                        color = "location") 

p$layers <- p$layers[-1]

# add aesthetics
p +
  geom_jitter(aes(shape = location), width = 0.2, size = 5, alpha = 0.8, show.legend = FALSE)+
  geom_vline(xintercept = 0, color="snow4") +
  geom_hline(yintercept = 0, color="snow4") +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons,
                     label = "p.signif", symnum.args = symnum.args) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  geom_boxplot(alpha = 0.7, col = "black") +
  scale_x_discrete(element_blank()) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18))


### determine whether the fractional contributions of glacial meltwater sourcing
### explains the variability in alpha diversity metrics 

# extract diversity metrics
palphadt <- data.table(p$data) 
palphadt_shan <- palphadt[(variable == "Shannon")]
palphadt_obs <- palphadt[(variable == "Observed")]


# Define the plot
a <- ggplot(data = palphadt_shan, 
            mapping = aes(fi_best, value,
                          color = location, shape = location)) +
  geom_point(size = 3, alpha = 0.8, show.legend = FALSE) + 
  facet_wrap(~location, ncol = 1) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  labs(y = "Shannon Index", x = "f<sub><i>i") +
  geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_markdown(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

b <- ggplot(data = palphadt_obs, 
            mapping = aes(fi_best, value,
                          color = location, shape = location)) +
  geom_point(size = 3, alpha = 0.8, show.legend = FALSE) + 
  facet_wrap(~location, ncol = 1) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  labs(y = "Observed Richness", x = "f<sub><i>i") +
  geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_markdown(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

plot_grid(b, a, labels = "auto", label_size = 20)




##### test comparisons if removing replicates $#####

# export alphadiversity files to manually average across 10 sites with replicates
write.csv(palphadt_shan, "shan_reps_test.csv")
write.csv(palphadt_obs, "obs_reps_test.csv")

# import an averaged dataset
shan_avg <- read.csv("shan_reps_test_avg.csv")
obs_avg <- read.csv("obs_reps_test_avg.csv")

# test the hypotheses that mean alpha diversity metrics differ between sites - For replicates
t.test(value ~ location, data = obs_avg, var.equal = TRUE) # for richness
wilcox.test(value ~ location, data = shan_avg) # for shannon diversity

# visualizing mean differences for metrics and sd differences
obs_avg %>% # for richness
  group_by(location) %>%
  summarize(sd = sd(value),
            mean = mean(value))

shan_avg %>% # for Shannon
  group_by(location) %>%
  summarize(sd = sd(value),
            mean = mean(value))
  
alphadt_shan_reps <- palphadt_shan %>% # selecting sites with replicates
  filter(sample_name == "Buzzard Fountain Spring" |
           sample_name == "Choss Seep" |
           sample_name == "Grinnell Spring" |
           sample_name == "James Spring" |
           sample_name == "Lunch Creek Spring" |
           sample_name == "Mineral Creek Tributary" |
           sample_name == "Mountain Goat Spring" |
           sample_name == "Palmer B Spring" |
           sample_name == "Paradise Park Spring" |
           sample_name == "Piegan North Spring" |
           sample_name == "Piegan Spring")

# visualize shannon diversity differences between sites with multiple replicates
alphadt_shan_reps %>%
  ggplot(aes(x = bio_rep, y = value)) +
  #geom_boxplot(outlier.shape  = NA) +
  geom_point(aes(color = sample_name), size = 6, alpha = 0.9, height = 0, width = .2) +
  #scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Shannon's Diversity", color = "Sample Site") +
  scale_y_continuous(limits = c(0,8)) +
  facet_wrap(~ sample_name) +
  theme_light()

# visualize OTU richness differences between sites with multiple replicates
alphadt_obs_reps <- palphadt_obs %>%
  filter(sample_name == "Buzzard Fountain Spring" |
           sample_name == "Choss Seep" |
           sample_name == "Grinnell Spring" |
           sample_name == "James Spring" |
           sample_name == "Lunch Creek Spring" |
           sample_name == "Mountain Goat Spring" |
           sample_name == "Mineral Creek Tributary" |
           sample_name == "Palmer B Spring" |
           sample_name == "Paradise Park Spring" |
           sample_name == "Piegan North Spring" |
           sample_name == "Piegan Spring")


alphadt_obs_reps %>%
  ggplot(aes(x = bio_rep, y = value)) +
  #geom_boxplot(outlier.shape  = NA) +
  geom_point(aes(color = sample_name), size = 6, alpha = 0.9, height = 0, width = .2) +
  #scale_color_manual(values = CPCOLS) +
  labs(x = "", y = "Observed", color = "Sample Site") +
  #scale_y_continuous(limits = c(0,8)) +
  facet_wrap(~ sample_name) +
  theme_light()

####### End replicate testing ######

## extract data for hypothesis testing
# separate by location and diversity metric
palphadt_shan_mh <- filter(palphadt_shan, location == "Mount Hood")
palphadt_shan_gnp <- filter(palphadt_shan, location == "Glacier National Park")
palphadt_obs_mh <- filter(palphadt_obs, location == "Mount Hood")
palphadt_obs_gnp <- filter(palphadt_obs, location == "Glacier National Park")

# fi vs. shannon's diversity index
# for MH
cor.test(x = palphadt_shan_mh %>% pull(fi_best),
         y = palphadt_shan_mh %>% pull(value)) %>%
  tidy()

# for GNP
cor.test(x = palphadt_shan_gnp %>% pull(fi_best),
         y = palphadt_shan_gnp %>% pull(value)) %>%
  tidy()

# fi vs. richness
# for MH
cor.test(x = palphadt_obs_mh %>% pull(fi_best),
         y = palphadt_obs_mh %>% pull(value)) %>%
  tidy()

# for GNP
cor.test(x = palphadt_obs_gnp %>% pull(fi_best),
         y = palphadt_obs_gnp %>% pull(value)) %>%
  tidy()

palphadt_obs %>%
  summarise(mh_range = )
browseVignettes("DESeq2")

# test whether mean fi differs between regions
t.test(fi_best ~ location, data = palphadt_obs)

# compare GNP and MH correlations via permutation -- Shannon
obs_cors_shan <- palphadt_shan %>%
  group_by(location) %>%
  summarise(cor = cor(fi_best,value)) %>%
  summarise(diff_cor = diff(cor)) %>%
  pull()

perm_shan <- replicate(10000, simplify = FALSE,
                       expr = palphadt_shan %>% 
                         mutate(perm_location = sample(location, replace = FALSE)) %>%
                         group_by(perm_location) %>%
                         summarise(cor = cor(fi_best,value)) %>%
                         summarise(diff_cor = diff(cor))) %>%
  bind_rows()

summarise(perm_shan, pval = mean(abs(diff_cor) >= abs(obs_cors_shan)))

# compare GNP and MH correlations via permutation -- OTU richness
obs_cors_obs <- palphadt_obs %>%
  group_by(location) %>%
  summarise(cor = cor(fi_best,value)) %>%
  summarise(diff_cor = diff(cor)) %>%
  pull()

perm_shan <- replicate(10000, simplify = FALSE,
                       expr = palphadt_obs %>% 
                         mutate(perm_location = sample(location, replace = FALSE)) %>%
                         group_by(perm_location) %>%
                         summarise(cor = cor(fi_best,value)) %>%
                         summarise(diff_cor = diff(cor))) %>%
  bind_rows()

summarise(perm_shan, pval = mean(abs(diff_cor) >= abs(obs_cors_shan)))
