# load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(phyloseq)
library(wesanderson)
library(vegan)
library(DESeq2)
library(gridExtra)
library(GGally)


# reference(s)
# Miller, J. B., Frisbee, M. D., Hamilton, T. L., & Murugapiran, S. K. (2021). 
# Recharge from glacial meltwater is critical for alpine springs and their microbiomes. 
# Environmental Research Letters, 16(6), 64012-. https://doi.org/10.1088/1748-9326/abf06b 

#code sourcing (partial): https://rpubs.com/lconteville/713954

# import mothur output produced from Miller et al. (2021)
sharedfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile <- "all18S_trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
mapfile <- "Master_Data_Sheet_11112021_consolidated_18s.csv"

# returns a phyloseq object with taxonomy and read counts where the taxa are rows
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
# import metadata
map <- read.csv(mapfile)

# rename location values to improve figure clarity
map$location <- recode(map$location, MH = "Mount Hood") %>% 
  recode(GNP = "Glacier National Park") 

# filter to extract samples of interest and avoid duplicate env values from replicates
map_filt <- filter(map, sample_type == "Spring" | sample_type == "Stream") %>%
  filter(ncbi_name_18s != "")
map_filt <- sample_data(map_filt)
rownames(map_filt) <- map_filt$ncbi_name_18s


# merge mothur data with metadata to create phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map_filt)

# replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# filter data (this is a choice depending on the project).
# below only takes otus with 1+ occurrences and <10% of total reads
# filter.data <- filter_taxa(data, function(x) sum(x > 0) > (0.1 * length(x)), TRUE)

# remove singletons
data_filt <- filter_taxa(moth_merge, function (x) {sum(x > 0) > 1}, prune=TRUE)

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
rare <- rarecurve(otu_tab, step=100, ylab="OTU",  label=F, col = cols, lty = line)

# how many reads in each sample?
row_sums <- rowSums(otu_tab)
sorted_row_sums <- sort(row_sums)
print(sorted_row_sums)


# identify which samples had fewer than 6000 reads and remove (12)
data.filt.clean <- subset_samples(data_filt, sample_sums(data_filt) >= (5000))
min(sample_sums(data.filt.clean)) #check min at 8500. revise above to be 90% of this min (8792)

# how many samples remain from each mountain range? --> threshold created to align with alpha-diversity sample size
sample_data(data.filt.clean) %>%
  group_by(location) %>%
  summarise(n = n())
df_data <- data.frame(sample_data(data.filt.clean))
write.csv(df_data$ncbi_name_16s, "sam_scope.csv")

## informative guide: http://joey711.github.io/phyloseq-extensions/DESeq2.html
vst.data.deseq <- phyloseq_to_deseq2(data.filt.clean, ~1)
# error:every gene contains at least one zero, cannot compute log geometric means
# fix: https://github.com/joey711/phyloseq/issues/387 & https://github.com/joey711/phyloseq/issues/445 , since several OTUs have 0's
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(vst.data.deseq), 1, gm_mean)
vst.data.deseq = estimateSizeFactors(vst.data.deseq, geoMeans=geoMeans)
vst.data.deseq <- estimateDispersions(vst.data.deseq, fitType = "local")

sum(otu_table(vst.data.deseq)) # n_reads
nsamples(data.filt.clean) # n_samples
ntaxa(data.filt.clean) # n_OTUs

#copy ps and add transformed counts, change negative counts to 0
vst.data <- data.filt.clean
otu_table(vst.data) <- otu_table(getVarianceStabilizedData(vst.data.deseq),
                                 taxa_are_rows = TRUE)
min(vst.data@otu_table) #-2.19

### because negative values after transformation likely
### represent zero counts or very few counts and for the
### distances and hypothesis in future tests these values would 
### be negligible or inconsequential

vst.data.0 <- transform_sample_counts(vst.data, function(x) {
  x[x < 0] <- 0
  x })

sum(otu_table(vst.data.0)) # n_reads
nsamples(vst.data.0) # n_samples
ntaxa(vst.data.0) # n_OTUs

sample_data(vst.data.0) %>%
  group_by(location) %>%
  summarise(n = n())

reps_explore <- sample_data(vst.data.0)

#How many counts?
sum(sample_sums(vst.data.0))

#ORDINATE - Regional
pca_all <- ordinate(vst.data.0, "RDA", correlation=TRUE)
#pull components to scale plot by top 2 causes of variance
pca_all_pc1 <- pca_all$CA$eig[1] / sum(pca_all$CA$eig)
pca_all_pc2 <- pca_all$CA$eig[2] / sum(pca_all$CA$eig)

#REPEAT for MH and GNP, exclusively. Will be used later on
vst_mh <- subset_samples(vst.data.0, location == "Mount Hood")
vst_mh_nowm <- subset_samples(vst.data.0, location == "Mount Hood"
                              & region != "Warm Springs")
vst_gnp <- subset_samples(vst.data.0, location =="Glacier National Park")

pca_mh <- ordinate(vst_mh, "RDA", correlation=TRUE)
pca_mh_nowm <- ordinate(vst_mh_nowm, "RDA", correlation=TRUE)
pca_gnp <- ordinate(vst_gnp, "RDA", correlation=TRUE)

pca_mh_pc1 <- pca_mh$CA$eig[1] / sum(pca_mh$CA$eig)
pca_mh_pc2 <- pca_mh$CA$eig[2] / sum(pca_mh$CA$eig)

pca_mh_nowm_pc1 <- pca_mh_nowm$CA$eig[1] / sum(pca_mh_nowm$CA$eig)
pca_mh_nowm_pc2 <- pca_mh_nowm$CA$eig[2] / sum(pca_mh_nowm$CA$eig)

pca_gnp_pc1 <- pca_gnp$CA$eig[1] / sum(pca_gnp$CA$eig)
pca_gnp_pc2 <- pca_gnp$CA$eig[2] / sum(pca_gnp$CA$eig)

plot_ordination(vst.data.0, pca_all, type="sites", color="location", 
                shape = "location") +
  #defaults from HS
  geom_vline(xintercept = 0, color="snow4") +
  geom_hline(yintercept = 0, color="snow4") +
  geom_point(size = 5, alpha=0.7) +
  #add ellipses. dotted = normal distribution; line = t-distribution
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") + 
  coord_fixed((pca_all_pc2) / pca_all_pc1) + 
  scale_color_manual(values=wes_palette("Moonrise2"),
                     labels=c("Glacier National Park", "Mount Hood")) +
  scale_shape_manual(values=c(16,17))+
  #scale_shape_manual(values = c(17,16), 
  #breaks ="location",
  #labels=c("Glacier National Park","Mount Hood"), 
  #name="Location") +
  theme_light() +
  theme(legend.position = "top")

plot_ordination(vst_gnp, pca_gnp, type="sites", color="region") +
  #defaults from HS
  geom_vline(xintercept = 0, color="snow4") +
  geom_hline(yintercept = 0, color="snow4") +
  geom_point(size = 5, alpha=0.7) +
  #add ellipses. dotted = normal distribution; line = t-distribution
  #stat_ellipse(type = "norm", linetype = 2) +
  #stat_ellipse(type = "t") + 
  coord_fixed((pca_gnp_pc2) / pca_gnp_pc1) + 
  #scale_color_manual(values=wes_palette("Moonrise2"),
                     #labels=c("Glacier National Park", "Mount Hood")) +
  #scale_shape_manual(values=c(16,17))+
  #scale_shape_manual(values = c(17,16), 
  #breaks ="location",
  #labels=c("Glacier National Park","Mount Hood"), 
  #name="Location") +
  theme_light() +
  theme(legend.position = "right")

#PERMANOVA
set.seed(1)

# Calculate bray curtis distance matrix
vst_bray <- phyloseq::distance(vst.data.0, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(vst.data.0))

# Adonis test
adonis2(vst_bray ~ location, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(vst_bray, sampledf$location)
permutest(beta)

#PERMANOVA - GNP
set.seed(1)

# Calculate bray curtis distance matrix
vst_bray_gnp <- phyloseq::distance(vst_gnp, method = "bray")

# make a data frame from the sample_data
sampledf_gnp <- data.frame(sample_data(vst_gnp))

# Adonis test
adonis(vst_bray_gnp ~ region, data = sampledf_gnp)

# Homogeneity of dispersion test
beta <- betadisper(vst_bray_gnp, sampledf_gnp$region)
permutest(beta)

#PERMANOVA - MH_NOWM

# Calculate bray curtis distance matrix
vst_bray_mh_nowm <- phyloseq::distance(vst_mh_nowm, method = "bray")

# make a data frame from the sample_data
sampledf_mh_nowm <- data.frame(sample_data(vst_mh_nowm))

# Adonis test
adonis(vst_bray_mh_nowm ~ region, data = sampledf_mh_nowm)

# Homogeneity of dispersion test
beta <- betadisper(vst_bray_mh_nowm, sampledf_mh_nowm$region)
permutest(beta)

#PERMANOVA - MH

# Calculate bray curtis distance matrix
vst_bray_mh <- phyloseq::distance(vst_mh, method = "bray")

# make a data frame from the sample_data
sampledf_mh <- data.frame(sample_data(vst_mh))

# Adonis test
adonis(vst_bray_mh ~ region, data = sampledf_mh)

# Homogeneity of dispersion test
beta <- betadisper(vst_bray_mh, sampledf_mh$region)
permutest(beta)

#Fitting Environmental Vectors

# Extract metadata. Need to be the same length as the PCA
# data, but doesn't relate to the taxonomy or shared files
# meta cheatsheet: fi_best, elevation_m, ysi_temp, ysi_do, ysi_spc,
#ysi_ph, ysi_orp, alkalinity_caco3, bicarbonate_hco3,
#calcium, chloride, magnesium, potassium, sio2, sodium,
#strontium, sulfate, tds
meta_table <- data.frame(fi = (sample_data(vst.data.0)$fi_best)^2,
                         Elevation = (sample_data(vst.data.0)$elevation_m),
                         Temp = log10(sample_data(vst.data.0)$ysi_temp),
                         DO = (sample_data(vst.data.0)$ysi_do)^2,
                         SPC = log10(sample_data(vst.data.0)$ysi_spc),
                         pH = (sample_data(vst.data.0)$ysi_ph),
                         ORP = (sample_data(vst.data.0)$ysi_orp),
                         CaCO3 = log10(sample_data(vst.data.0)$alkalinity_caco3),
                         HCO3 = log10(sample_data(vst.data.0)$bicarbonate_hco3),
                         Ca = log10(sample_data(vst.data.0)$calcium),
                         Cl = log10(sample_data(vst.data.0)$chloride),
                         Mg = log10(sample_data(vst.data.0)$magnesium),
                         K = log10(sample_data(vst.data.0)$potassium),
                         SiO2 = log10(sample_data(vst.data.0)$sio2),
                         Na = log10(sample_data(vst.data.0)$sodium),
                         Sr = log10(sample_data(vst.data.0)$strontium),
                         SO4 = log10(sample_data(vst.data.0)$sulfate),
                         TDS = log10(sample_data(vst.data.0)$tds))

## Check correlations and skew. Then transform variables in 
## meta.table as needed
#left = scatterplots between two different variables. right = pearson correlation.
#middle = variable distribution
#ggpairs(meta_table)
###how do i know what to fix, given the result of this gg pairs figure?###

# Fit Environmental metadata to 16s PCA
set.seed(503)
env_fit <- envfit(pca_all, env = meta_table,
                       na.rm = TRUE, permutations = 9999, display = "sites")
env_fit_adj <- env_fit 
env_fit_adj_pvals <- p.adjust(env_fit$vectors$pvals, method = "bonferroni") 
env_fit_adj$vectors$pvals <- env_fit_adj_pvals
env_fit_adj

library(factoextra)
library(ellipse)

palette(c(adjustcolor("#667c74", alpha.f=0.6), adjustcolor("#b3692b", alpha.f = 0.6)))
par(mar=c(5,6,4,1)+.1)
p <- ordiplot(pca_all, display = "sites", type = "n", xlab = "PC1 [17.4%]", ylab = "PC2 [10.8%]", #find pc % explained by plotting in ggplot
              xlim=c(-10,20), ylim=c(-10,10), cex.lab = 2)
plot(env_fit_adj, p = 0.05, add = TRUE, col = "slategrey")
points(pca_all, display = "sites", pch = c(16,17)[as.factor(vst.data.0@sam_data$location)],
       col = as.factor(vst.data.0@sam_data$location), cex = 2)
ordiellipse(p, vst.data.0@sam_data$location, kind = "sd", conf=0.95, 
            col = c(adjustcolor("#667c74", alpha.f=0.6), adjustcolor("#b3692b", alpha.f = 0.6)),
            lwd=3)
#ordiellipse(p, vst.data.0@sam_data$location, kind = "ehull", conf=0.95, 
            #col = c(adjustcolor("#667c74", alpha.f=0.6), adjustcolor("#b3692b", alpha.f = 0.6)),
            #draw = "polygon", lty = 0)
legend(x = 9, y = 10, inset = 0.02, col=unique(as.factor(vst.data.0@sam_data$location)),
                                                          pch= c(17, 16), 
       legend = unique(as.factor(vst.data.0@sam_data$location)), box.lty = 0, cex = 1.3)

#REPEAT fFitting Var for MH
meta_table_mh <- data.frame(fi = (sample_data(vst_mh)$fi_best),
                         Elevation = (sample_data(vst_mh)$elevation_m),
                         Temp = log10(sample_data(vst_mh)$ysi_temp),
                         DO = (sample_data(vst_mh)$ysi_do)^2,
                         #SPC = log10(sample_data(vst_mh)$ysi_spc), -> redundant w/ mg group
                         #pH = (sample_data(vst_mh)$ysi_ph), -> redundant w/ mg group
                         ORP = log10(sample_data(vst_mh)$ysi_orp),
                         #CaCO3 = log10(sample_data(vst_mh)$alkalinity_caco3), -> redundant w/ mg group
                         #HCO3 = log10(sample_data(vst_mh)$bicarbonate_hco3), -> redundant w/ mg group
                         #Ca = log10(sample_data(vst_mh)$calcium), -> redundant w/ mg group
                         #Cl = log10(sample_data(vst_mh)$chloride),-> redundant w/ T group
                         Mg = log10(sample_data(vst_mh)$magnesium),
                         K = log10(sample_data(vst_mh)$potassium),
                         SiO2 = log10(sample_data(vst_mh)$sio2),
                         Na = log10(sample_data(vst_mh)$sodium),
                         Sr = log10(sample_data(vst_mh)$strontium),
                         SO4 = log10(sample_data(vst_mh)$sulfate),
                         #TDS = log10(sample_data(vst_mh)$tds), -> redundant w/ Na
                         Dispersal = log10(sample_data(vst_mh)$dist_from_summit))

meta_table_mh_nowm <- data.frame(fi = (sample_data(vst_mh_nowm)$fi_best),
                                 Elevation = (sample_data(vst_mh_nowm)$elevation_m),
                                 Temp = log10(sample_data(vst_mh_nowm)$ysi_temp),
                                 DO = (sample_data(vst_mh_nowm)$ysi_do)^2,
                                 SPC = log10(sample_data(vst_mh_nowm)$ysi_spc),
                                 pH = (sample_data(vst_mh_nowm)$ysi_ph),
                                 ORP = log10(sample_data(vst_mh_nowm)$ysi_orp),
                                 #CaCO3 = log10(sample_data(vst_mh_nowm)$alkalinity_caco3), -> redundant w/ hco3
                                 HCO3 = log10(sample_data(vst_mh_nowm)$bicarbonate_hco3),
                                 Ca = log10(sample_data(vst_mh_nowm)$calcium),
                                 Cl = log10(sample_data(vst_mh_nowm)$chloride),
                                 Mg = log10(sample_data(vst_mh_nowm)$magnesium),
                                 K = log10(sample_data(vst_mh_nowm)$potassium),
                                 SiO2 = log10(sample_data(vst_mh_nowm)$sio2),
                                 Na = log10(sample_data(vst_mh_nowm)$sodium),
                                 Sr = log10(sample_data(vst_mh_nowm)$strontium),
                                 SO4 = log10(sample_data(vst_mh_nowm)$sulfate),
                                 TDS = log10(sample_data(vst_mh_nowm)$tds),
                                 Dispersal = log10(sample_data(vst_mh_nowm)$dist_from_summit))

## Check correlations and skew. Then transform variables in 
## meta.table as needed
#left = scatterplots between two different variables. right = pearson correlation.
#middle = variable distribution
ggpairs(meta_table_mh)
ggpairs(meta_table_mh_nowm)
###how do i know what to fix, given the result of this gg pairs figure?###

# Fit Environmental metadata to 16s PCA
set.seed(503)
env_fit_mh <- envfit(pca_mh, env = meta_table_mh,
                  na.rm = TRUE, permutations = 9999, display = "sites")
env_fit_adj_mh <- env_fit_mh
env_fit_adj_pvals_mh <- p.adjust(env_fit_mh$vectors$pvals, method = "bonferroni") 
env_fit_adj_mh$vectors$pvals <- env_fit_adj_pvals_mh
env_fit_adj_mh

#no wm
set.seed(503)
env_fit_mh_nowm <- envfit(pca_mh_nowm, env = meta_table_mh_nowm,
                     na.rm = TRUE, permutations = 9999, display = "sites")
env_fit_adj_mh_nowm <- env_fit_mh_nowm
env_fit_adj_pvals_mh_nowm <- p.adjust(env_fit_mh_nowm$vectors$pvals, method = "bonferroni") 
env_fit_adj_mh_nowm$vectors$pvals <- env_fit_adj_pvals_mh_nowm
env_fit_adj_mh_nowm

#need to know variance explained by PCA
plot_ordination(vst_mh, pca_mh, type="sites", color="region", 
                shape = "region")

#w/o wm
plot_ordination(vst_mh_nowm, pca_mh_nowm, type="sites", color="region", 
                shape = "region")

library(factoextra)
library(ellipse)

palette(c(adjustcolor("#FE938C", alpha.f=0.6), 
          adjustcolor("#E6B89C", alpha.f = 0.6),
          adjustcolor("#4281A4", alpha.f = 0.6),
          adjustcolor("#50514F", alpha.f = 0.6),
          adjustcolor("#ff1654", alpha.f = 0.6),
          adjustcolor("#70C1B3", alpha.f = 0.6)))
          
p_mh <- ordiplot(pca_mh, display = "sites", type = "n", xlab = "PC1 [13.4%]", ylab = "PC2 [12.3%]",
             xlim = c(-10,20), ylim = c(-12,15),
              cex.lab = 2)
plot(env_fit_adj_mh, p = 0.05, add = TRUE, col = "slategrey")
points(pca_mh, display = "sites",
       col = as.factor(vst_mh@sam_data$region), pch = 16, cex = 2)
ordiellipse(p_mh, vst_mh@sam_data$region, kind = "sd", conf=0.95, 
            col = c(adjustcolor("#FE938C", alpha.f=0.6), 
                    adjustcolor("#E6B89C", alpha.f = 0.6),
                    adjustcolor("#4281A4", alpha.f = 0.6),
                    adjustcolor("#50514F", alpha.f = 0.6),
                    adjustcolor("#ff1654", alpha.f = 0.6),
                    adjustcolor("#70C1B3", alpha.f = 0.6)),
            lwd=3)
legend(x = 12.5, y = 16, inset = 0.02, pch = 16, col=unique(as.factor(vst_mh@sam_data$region)), 
       legend = unique(as.factor(vst_mh@sam_data$region)), box.lty = 0, cex = 1.3)

#w/o wm
palette(c(adjustcolor("#FE938C", alpha.f=0.6), 
          adjustcolor("#E6B89C", alpha.f = 0.6),
          adjustcolor("#4281A4", alpha.f = 0.6),
          adjustcolor("#50514F", alpha.f = 0.6),
          adjustcolor("#70C1B3", alpha.f = 0.6)))

p_mh_nowm <- ordiplot(pca_mh_nowm, display = "sites", type = "n", 
                      xlab = "PC1 [15.4%]", ylab = "PC2 [12%]",
                      ylim = c(-12,14), xlim = c(-15, 20), cex.lab = 2)
plot(env_fit_adj_mh_nowm, p = 0.05, add = TRUE, col = "slategrey")
points(pca_mh_nowm, display = "sites",
       col = as.factor(vst_mh_nowm@sam_data$region), pch = 16, cex = 2)
ordiellipse(p_mh_nowm, vst_mh_nowm@sam_data$region, kind = "sd", conf=0.95, 
            col = c(adjustcolor("#FE938C", alpha.f=0.6), 
                    adjustcolor("#E6B89C", alpha.f = 0.6),
                    adjustcolor("#4281A4", alpha.f = 0.6),
                    adjustcolor("#50514F", alpha.f = 0.6),
                    adjustcolor("#70C1B3", alpha.f = 0.6)),
            lwd=3)
legend(x = 14, y = 14, inset = 0.02, pch = 16, col=unique(as.factor(vst_mh_nowm@sam_data$region)), 
       legend = unique(as.factor(vst_mh_nowm@sam_data$region)), box.lty = 0, cex = 1.3)


#REPEAT Fitting Var for GNP
meta_table_gnp <- data.frame(Fi = (sample_data(vst_gnp)$fi_best),
                            Elevation = (sample_data(vst_gnp)$elevation_m)^2,
                            Temp = log10(sample_data(vst_gnp)$ysi_temp),
                            DO = (sample_data(vst_gnp)$ysi_do),
                            SPC = (sample_data(vst_gnp)$ysi_spc),
                            pH = (sample_data(vst_gnp)$ysi_ph),
                            ORP = (sample_data(vst_gnp)$ysi_orp),
                            CaCO3 = log10(sample_data(vst_gnp)$alkalinity_caco3),
                            HCO3 = (sample_data(vst_gnp)$bicarbonate_hco3),
                            Ca = (sample_data(vst_gnp)$calcium),
                            Cl = log10(sample_data(vst_gnp)$chloride),
                            Mg = log10(sample_data(vst_gnp)$magnesium),
                            K = (sample_data(vst_gnp)$potassium),
                            SiO2 = (sample_data(vst_gnp)$sio2),
                            Na = (sample_data(vst_gnp)$sodium),
                            Sr = (sample_data(vst_gnp)$strontium),
                            SO4 = (sample_data(vst_gnp)$sulfate),
                            TDS = (sample_data(vst_gnp)$tds))

## Check correlations and skew. Then transform variables in 
## meta.table as needed
#left = scatterplots between two different variables. right = pearson correlation.
#middle = variable distribution
ggpairs(meta_table_gnp)
###how do i know what to fix, given the result of this gg pairs figure?###

# Fit Environmental metadata to 16s PCA
set.seed(503)
env_fit_gnp <- envfit(pca_gnp, env = meta_table_gnp,
                     na.rm = TRUE, permutations = 9999, display = "sites")
env_fit_adj_gnp <- env_fit_gnp
env_fit_adj_pvals_gnp <- p.adjust(env_fit_gnp$vectors$pvals, method = "bonferroni") 
env_fit_adj_gnp$vectors$pvals <- env_fit_adj_pvals_gnp
env_fit_adj_gnp

#need to know variance explained by PCA
plot_ordination(vst_gnp, pca_gnp, type="sites", color="region", 
                shape = "region")

library(factoextra)
library(ellipse)

palette(c(adjustcolor("#FE938C", alpha.f=0.6), 
          adjustcolor("#EAD2AC", alpha.f = 0.6),
          adjustcolor("#9CAFB7", alpha.f = 0.6),
          adjustcolor("#4281A4", alpha.f = 0.6),
          adjustcolor("#50514F", alpha.f = 0.6),
          adjustcolor("#ff1654", alpha.f = 0.6),
          adjustcolor("#70C1B3", alpha.f = 0.6)))

p_gnp <- ordiplot(pca_gnp, display = "sites", type = "n", xlab = "PC1 [21.9%]", ylab = "PC2 [14.2%]", 
                  cex.lab = 2)
                 #xlim = c(-5,5), ylim = c(-10,10))
plot(env_fit_adj_gnp, p = 0.05, add = TRUE, col = "slategrey")
points(pca_gnp, display = "sites",
       col = as.factor(vst_gnp@sam_data$region), pch = 16, cex = 3)
#ordiellipse(p_gnp, vst_gnp@sam_data$region, kind = "sd", conf=0.95, 
            #col = c(adjustcolor("#FE938C", alpha.f=0.6), 
                    #adjustcolor("#EAD2AC", alpha.f = 0.6),
                    #adjustcolor("#9CAFB7", alpha.f = 0.6),
                    #adjustcolor("#4281A4", alpha.f = 0.6),
                    #adjustcolor("#50514F", alpha.f = 0.6),
                    #adjustcolor("#ff1654", alpha.f = 0.6),
                    #adjustcolor("#70C1B3", alpha.f = 0.6)),
            #lwd=3)
legend(x = -10, y = 2, pch = 16, inset = 0.02, col=unique(as.factor(vst_gnp@sam_data$region)), 
       legend = unique(as.factor(vst_gnp@sam_data$region)), box.lty = 0, cex = 1.3)

