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

#### Phyloseq Object Setup) ####

#Import data from MSI: /panfs/roc/groups/6/hamil689/shared/sen/JM_second_paper/res/16SII
sharedfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated_18s.csv"


#Returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

colnames <- colnames(otu_table(mothur_data))

#Import metadata -> remove blanks and samples lacking recharge_source
map <- read.csv(mapfile)
map$location <- recode(map$location, MH = "Mount Hood") %>% #rename location values to improve figure clarity
  recode(GNP = "Glacier National Park")


#visual spread of fi best
ggplot(data = map, mapping = aes(x = location, y = fi_best, 
                                 color = location, shape = location)) +
  geom_jitter(size = 5, alpha = 0.6, width = 0.2, show.legend = FALSE) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  theme_light()


#remove blanks
map.recharge <- filter(map, bio_rep != "")
map.recharge <- filter(map.recharge, ncbi_name_18s != "")
map.recharge <- mutate(map.recharge, primary_source = 
                         case_when(fi_best > 0.6 ~ "Glacier Meltwater",
                                   fi_best <= 0.6 ~ "Precipitation"))
map.recharge <- sample_data(map.recharge)
rownames(map.recharge) <- map.recharge$ncbi_name_18s

#merge mothur data with metadata returns returns to create phyloseq object with otus, taxa and metadata are rows
moth_merge <- merge_phyloseq(mothur_data, map.recharge)

sum(otu_table(moth_merge)) # n_reads
nsamples(moth_merge) # n_samples
ntaxa(moth_merge) # n_OTUs



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

#provides alpha diversity estimates for each sample
#estimate_richness(data_filt_clean_srs)
a <- estimate_richness(data_filt_rarified)
sum(sample_sums(data_filt_rarified))
    
comparisons <- list( c("Glacier National Park", "Mount Hood"))
args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")) #manually add wilcox values dervied below for aesthetics

#plot results from above estimated alpha diversity. isolating observed and shannon
p <- plot_richness(data_filt_rarified, x="location", 
              measures=c("Observed", "Shannon"),
              color = "location") 

p$layers <- p$layers[-1]
p +
  geom_jitter(aes(shape = location), width = 0.2, size = 5, alpha = 0.6, show.legend = FALSE)+
  geom_vline(xintercept = 0, color="snow4") +
  geom_hline(yintercept = 0, color="snow4") +
  stat_compare_means(method = "t.test", comparisons = comparisons,
                     label = "p.signif", symnum.args = args) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  stat_summary(fun.data = "mean_cl_boot",
               show.legend = FALSE,
               position = position_nudge(x = .3, y = 0),
               fun.args=(conf.int=0.95),
               color = "black",
               size = 0.7) +
  scale_x_discrete(element_blank()) +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        strip.text.x = element_text(size = 18))

#To determine correct two-sample test, review whether assumptions are met for both Shannon and Observed
#Extract diversity metrics by sample
alpha_vals <- select(p$data, location, sample_name, variable, value)
alpha_shannon <- filter(alpha_vals, variable == "Shannon")
alpha_richness <- filter(alpha_vals, variable == "Observed")
#Shannon vars equal? -> Yes! Within factor of 5
alpha_shannon %>%
  group_by(location) %>%
  summarise(mean = mean(value),
            var = var(value),
            sd = sd(value))
#Richness vars equal? -> Yes! Within factor of 5
alpha_richness %>%
  group_by(location) %>%
  summarise(mean = mean(value),
            var = var(value),
            sd = sd(value))

#### Shannon vars normal? No per Q-Q plot. Therefore, we must permute
alpha_shannon_lm <- lm(value ~ location, data = alpha_shannon)
autoplot(alpha_shannon_lm)

#Conduct permutation
set.seed(503)

#Collect observed mean for both locations and generate observed mean-diff test statistic
obs <- alpha_shannon %>%
  group_by(location) %>%
  summarise(mean = mean(value)) %>%
  summarise(mean_diff = diff(mean)) %>%
  pull()

#Permute to generate a null distribution of mean differences
perm <- replicate(9999, simplify = FALSE,
                     expr = alpha_shannon %>%
                       mutate(perm_location = sample(location, replace = FALSE)) %>%
                       group_by(perm_location) %>%
                       summarise(mean = mean(value)) %>%
                       summarise(mean_diff = diff(mean))) %>%
  bind_rows()

#Derive p-value representing the probability of median-difference being more extreme than our observation
#under the null distribution
summarise(perm, pval = mean(abs(mean_diff) >= abs(obs)))

#Does t-test differ? Very minimally. Let's stick with it so we can use Cohens
wilcox.test(value ~ location, data = alpha_shannon, var.equal = TRUE)
#0.0346
shan_lm <- lm(value ~ location, data = alpha_shannon)
summary.lm(shan_lm)

#Determine richness uncertainty, but we'll need summary statistics first
alpha_shannon_summaries <- alpha_shannon %>%
  group_by(location) %>%
  summarise(n         = n(),
            df        = n - 1,
            mean = mean(value),
            var  = var(value))
alpha_shannon_summaries %>%
  summarise(est_diff        = abs(diff(mean)),
            pooled_variance = sum(df * var) / sum(df),
            se              = sqrt(sum(pooled_variance/n)),
            crit_t_95       = qt(p = 0.05/2, df = sum(df), lower.tail = FALSE),
            lower_95CI       = est_diff  - se * crit_t_95,
            upper_95CI       = est_diff  + se * crit_t_95)

#### Richness vars normal? Yes, therefore we can performa  t-test
alpha_richness_lm <- lm(value ~ location, data = alpha_richness)
autoplot(alpha_richness_lm)
t.test(value ~ location, data = alpha_richness, var.equal = TRUE)
#0.0004972#
rich_lm <- lm(value ~ location, data = alpha_richness)
summary.lm(rich_lm)
#Determine richness uncertainty, but we'll need summary statistics first
alpha_richness_summaries <- alpha_richness %>%
  group_by(location) %>%
  summarise(n         = n(),
            df        = n - 1,
            mean = mean(value),
            var  = var(value))
alpha_richness_summaries %>%
  summarise(est_diff        = abs(diff(mean)),
            pooled_variance = sum(df * var) / sum(df),
            se              = sqrt(sum(pooled_variance/n)),
            crit_t_95       = qt(p = 0.05/2, df = sum(df), lower.tail = FALSE),
            lower_95CI       = est_diff  - se * crit_t_95,
            upper_95CI       = est_diff  + se * crit_t_95)


#create alpha diversity against continuous variable
palphadt <- data.table(p$data)
palphadt_shan <- palphadt[(variable == "Shannon")]
palphadt_obs <- palphadt[(variable == "Observed")]

# Order by Days - NOT NECESSARY(?)
palphadt_shan <- palphadt_shan[order(fi_best)][(is.na(se))]
palphadt_obs <- palphadt_obs[order(fi_best)][(is.na(se))]

library(data.table)
library(cowplot)

# Define the plot
 a <- ggplot(data = palphadt_shan, 
       mapping = aes(fi_best, value,
                     color = location, shape = location)) +
  # shape = ReportedAntibioticUsage)) + 
  geom_point(size = 5, alpha = 0.6, show.legend = FALSE) + 
 # geom_path() +
  #geom_point(data = alphadt[(ReportedAntibioticUsage == "Yes")], 
            # size = 8, alpha = 0.35) +
  facet_wrap(~location, ncol = 1) +
  scale_color_manual(values=wes_palette("Moonrise2")) +
  labs(y = "Shannon Index", x = bquote(~f[i])) +
  geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
  theme_light() +
   theme(legend.position = "none",
         axis.title = element_text(size = 20),
         strip.text.x = element_text(size = 18))

 b <- ggplot(data = palphadt_obs, 
             mapping = aes(fi_best, value,
                           color = location, shape = location)) +
   # shape = ReportedAntibioticUsage)) + 
   geom_point(size = 5, alpha = 0.6, show.legend = FALSE) + 
   # geom_path() +
   #geom_point(data = alphadt[(ReportedAntibioticUsage == "Yes")], 
   # size = 8, alpha = 0.35) +
   facet_wrap(~location, ncol = 1) +
   scale_color_manual(values=wes_palette("Moonrise2")) +
   labs(y = "Observed Richness", x = bquote(~f[i])) +
   geom_smooth(method = "lm", formula = y ~ x, se = T, show.legend = FALSE, color = "cornsilk4") +
   theme_light() +
   theme(legend.position = "none",
         axis.title = element_text(size = 20),
         strip.text.x = element_text(size = 18))
 
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
 
 #summarise correlation data. manually alter for either alpha diverstiy statistic
 #FIRST separate by location
 palphadt_shan_mh <- filter(palphadt_shan, location == "Mount Hood")
 palphadt_shan_gnp <- filter(palphadt_shan, location == "Glacier National Park")
 palphadt_obs_mh <- filter(palphadt_obs, location == "Mount Hood")
 palphadt_obs_gnp <- filter(palphadt_obs, location == "Glacier National Park")
   
 library(broom)
 
 #t-test for correlations - shan
 cor.test(x = palphadt_shan_mh %>% pull(fi_best),
          y = palphadt_shan_mh %>% pull(value)) %>%
   tidy()
 
 cor.test(x = palphadt_shan_gnp %>% pull(fi_best),
          y = palphadt_shan_gnp %>% pull(value)) %>%
   tidy()
 
 #t-test for correlations - shan
 cor.test(x = palphadt_obs_mh %>% pull(fi_best),
          y = palphadt_obs_mh %>% pull(value)) %>%
   tidy()
 
 cor.test(x = palphadt_obs_gnp %>% pull(fi_best),
          y = palphadt_obs_gnp %>% pull(value)) %>%
   tidy()
 
 #Compare GNP and MH correlations via permutation

#Use polynomial regression?
lm_shan <- lm(value ~ location*fi_best, palphadt_shan)
summary.lm(lm_shan)
 
 #Begin with SHANNON
 obs_cors_shan <- palphadt_shan %>%
   group_by(location) %>%
   summarise(cor = cor(fi_best,value)) %>%
   summarise(diff_cor = diff(cor)) %>%
   pull()
 
 perm_shan <- replicate(9999, simplify = FALSE,
                      expr = palphadt_shan %>% 
                        mutate(perm_location = sample(location, replace = FALSE)) %>%
                        group_by(perm_location) %>%
                        summarise(cor = cor(fi_best,value)) %>%
                        summarise(diff_cor = diff(cor))) %>%
   bind_rows()
 
 summarise(perm_shan, pval = mean(abs(diff_cor) >= abs(obs_cors_shan)))
 
 #Now OBSERVED
 obs_cors_obs <- palphadt_obs %>%
   group_by(location) %>%
   summarise(cor = cor(fi_best,value)) %>%
   summarise(diff_cor = diff(cor)) %>%
   pull()
 
 perm_obs <- replicate(9999, simplify = FALSE,
                        expr = palphadt_obs %>% 
                          mutate(perm_location = sample(location, replace = FALSE)) %>%
                          group_by(perm_location) %>%
                          summarise(cor = cor(fi_best,value)) %>%
                          summarise(diff_cor = diff(cor))) %>%
   bind_rows()
 
 summarise(perm_obs, pval = mean(abs(diff_cor) >= abs(obs_cors_obs)))
 
 #How different is fi_best between regions?
 t.test(fi_best ~ location, data = palphadt_shan)
citation("tidyverse")
?cor.test

mh <- subset_samples(data.filt.clean, location == "Mount Hood")
gnp <- subset_samples(data.filt.clean, location == "Glacier National Park")

sample_data(mh)$fi_best

##### LEFse Analysis #####
library(microbiomeMarker)

# create location-specific datasets
data_filt_rarified_mh <- subset_samples(data_filt_rarified, location == "Mount Hood")
data_filt_rarefied_gnp <- subset_samples(data_filt_rarified, location == "Glacier National Park")

tax <- tax_table(data_filt_rarified_mh)
colnames(tax) <- c("Kingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")  # if OTU is species-level
tax_table(data_filt_rarified_mh) <- tax
head(sample_data(data_filt_rarified_mh))

tax <- tax_table(data_filt_rarified)
colnames(tax) <- c("Kingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")  # if OTU is species-level
tax_table(data_filt_rarified) <- tax
head(sample_data(data_filt_rarified))



head(tax_table(data_filt_rarified_mh)[,"Species"], useNA = "ifany")



# mh sourcing differences
lef_out<-run_lefse(data_filt_rarified_mh, group = "primary_source", norm = "none", 
                   kw_cutoff = 0.01, lda_cutoff = 2.5)
lef_out
plot_ef_bar(lef_out)

rank_names(data_filt_rarified_mh)

table(sample_data(data_filt_rarified)$primary_source)


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

lef_out_reg <-run_lefse(data_filt_rarified_mh_nonw, group = "quadrant", norm = "none", 
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
