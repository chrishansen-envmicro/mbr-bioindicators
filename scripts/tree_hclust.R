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

# filter to extract samples of interest and avoid duplicate env values from replicates
map_filt <- filter(map, sample_type == "Spring" | sample_type == "Stream") %>%
  filter(ncbi_name_16s != "")

map_filt_mut <- mutate(map_filt, primary_source = 
                         case_when(fi_best > 0.6 ~ "Glacier Meltwater",
                                   fi_best <= 0.6 ~ "Precipitation"))
map_filt_mut <- sample_data(map_filt_mut)
rownames(map_filt_mut) <- map_filt_mut$ncbi_name_16s


# merge mothur data with metadata to create phyloseq object
moth_merge <- merge_phyloseq(mothur_data, map_filt_mut)

# replace rank with taxonomy
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# remove non-bacterial mitochondria and chloropast seq. should be removed during amp processing (mother), but this is a double-check
data <- moth_merge %>% subset_taxa(Family != "mitochondria" & Class != "Chloroplast")

# filter data (this is a choice depending on the project).
# below only takes otus with 1+ occurrences and <10% of total reads
filter.data <- filter_taxa(data, function(x) sum(x > 0) > (0.1 * length(x)), TRUE)

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

otu_tab <- otu_table(filter.data) # extract otu table
class(otu_tab) <- "matrix" # as.matrix() will do nothing. you get a warning here, but this is what we need to have
otu_tab <- t(otu_tab) # transpose observations to rows
rare <- rarecurve(otu_tab, step=100, ylab="OTU",  label=F, col = cols, lty = line)
abline(v=6000, lw = 1) # display selected cutoff across samples that balance 

# identify which samples had fewer than 6000 reads and remove (12)
data.filt.clean <- subset_samples(filter.data, sample_sums(filter.data) >= (6000))
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

#copy ps and add transformed counts, change negative counts to 0
vst.data <- filter.data
otu_table(vst.data) <- otu_table(getVarianceStabilizedData(vst.data.deseq),
                                 taxa_are_rows = TRUE)
min(vst.data@otu_table) #-2.19

###because negative values after transformation likely
### represent zero counts or very few counts and for the
### distances and hypothesis in future tests these values would 
### be negligible or inconsequential

vst.data.0 <- transform_sample_counts(vst.data, function(x) {
  x[x < 0] <- 0
  x })

#use bray-certis dissimilarity to summarize distances between ecological samples
#note: unweighted = presence/absence of taxa
#weighted: relative abundance of taxa + phylogenetic distance
#source: https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html

#create hclust
library(dendextend)
library(ggdendro)
library(data.table)

# Bray-Curtis dissimilarity instead of UniFrac
hclust.data <- phyloseq::distance(vst.data.0, method = "bray")  # data.16s is your phyloseq object

# Create df of sample metadata
vst_dat <- data.frame(sample_data(vst.data.0))

# Create hclust using Ward.D2 method
hclust <- hclust(hclust.data, method = "ward.D2")

# Convert to dendrogram
hclustd <- as.dendrogram(hclust)

#ease digestibility: http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
nodePar <- list(lab.cex = 0.5, pch = c(NA, 19), 
                cex = 1.8, col = "blue")

#need to change rownames to hydro_names + setup at tibble to dictate color, etc
dend_data <- dendro_data(hclustd, type = "rectangle")
labs <- dend_data$labels
labstib <- tibble(labs$label)
labstib <- dplyr::rename(labstib, "ncbi_name_16s" = "labs$label")
mapssel <- select(vst_dat,ncbi_name_16s, primary_source, alt_sample_name, sample_type, location)
labsfin <- left_join(labstib, mapssel, by = "ncbi_name_16s")

#change label names (setup) - hydro_names
vec = pull(labsfin, alt_sample_name)
newclustd = hclustd
labels(newclustd) <- vec
#text colors (setup) - sample_type
#type = select(labsfin, sample_type)
#type = mutate(type, sample_type = recode(.x = sample_type,
                                        # "Snow Algae" = "grey", "Snow" = "grey", "Glacial Ice" = "grey",
                                        # "Stream" = "black", "Spring" = "black", "Glacial Melt" = "grey"))
#typed = pull(type, sample_type)
#fill colors (setup) - MH vs. GNP
loc = select(labsfin, location)
loc = mutate(loc, location = recode(.x = location, "Mount Hood" = "#c27d38", "Glacier National Park" = "#798e87"))
locd = pull(loc, location)
#change shape (setup) - glacier (19) vs. rain/snow recharge (1)
rech = select(labsfin, primary_source)
rech = mutate(rech, primary_source = recode(.x = primary_source, 
                                                 "Glacier Meltwater" = 19, "Precipitation" = 1))
rechd = pull(rech, primary_source)
#manual cluster bars
clus1 = vec[1:12]
clus2 = vec[13:35]
clus3 = vec[36:48]
clus4 = vec[49:60]

pars = par()
par(mar = c(11, 3, 2, 2))
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
newclustd %>% set("labels_cex", 1) %>% 
  set("leaves_pch", rechd) %>% 
  set("leaves_cex", 2.2) %>% 
  set("leaves_col", locd) %>% 
  #set("labels_col", typed) %>% 
  set("branches_lwd", 1) %>%
  set("branches_lty", 1) %>%
  #raise.dendrogram(-.25) %>% 
  set("branches_col", "darkgray") %>% 
  plot()