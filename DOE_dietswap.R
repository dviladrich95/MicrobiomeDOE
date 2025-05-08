# importing libraries
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)

#upload data dietswaps from package phyloseq (data from mucosa sample analysis regarding different dietary treatments)
data("dietswap")
dietswap

# show the first ten rows (species) and 7 first columns (samples) of the otu.table
dietswap@otu_table@.Data[1:10,1:7]

# get the metadata (= variables associated to each sample, e.g. sampleID, diet, collection time point, sex, nationality etc.) from each samples

metadata <- as.data.frame(sample_data(dietswap))

# correction line, in order to tell R that the columns timepoint and timepoint.within.group are categorical variables and not numerical
metadata <- metadata %>%
mutate(
timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group)
)

# summary of your data  
summary(metadata)

# psmelt() converts a phyloseq object into a flat data frame format — by "melting" it — so it's easier to use with ggplot2 (because phyloseq data is hierarchical, thus doesn’t allow for easy handling with ggplot2)
diet <- psmelt(dietswap)

# visualize data, e.g. the percentage of Afroamericans (AFR) or American (AAM) from the sample being obese, overweight or lean
plot_frequencies(x=sample_data(dietswap),
Groups="bmi_group", Factor="nationality")

# filter out non-prevalent taxa from our samples (allow to only compare those taxa which are frequent within all samples)
core_diet <- core(x=dietswap,
detection=50, # threshold for presence of each OTU in a sample
prevalence=50/100) # threshold for the fraction of OTUs that exceeds detection threshold => means that a taxon must appear in at least 50% of the samples to be included. "Appear" means its abundance is ≥ the detection threshold (in this case, 50).

# logartimithic transformation of our data allows for normalization (in microbiome data a few taxa are highly abundant while others are rare, high variance). For some packages other forms of normalization apply, so not always useful (e.g. DESeq2)
diet_log <- transform_sample_counts(core_diet, function(x) log(1+x))

# calculate alpha diversity metrics
rich <- estimate_richness(core_diet)
rich

# you could also select only certain indicators, e.g. Shannon and Observed gives the count of taxa or OTUs in each sample
rich <- estimate_richness(core_diet, measures = c("Observed", "Shannon"))
rich

# plot the Shannon metrics according to bmi group categories. Plot richness already allows for calculating and plotting diversity indices. 
plot_richness(core_diet, x="bmi_group", measures= c("Shannon")) +
geom_boxplot() +
theme_classic() +
theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

#  Assess significance of differences in microbial diversity across bmi groups through statistical test
## convert richness data into a data frame for easier handling 
rich_df <- as.data.frame(rich)

# add the bmi_group to the richness data 
rich_df$bmi_group <- sample_data(core_diet)$bmi_group

# use Wilcoxon rank-sum test for observed richness 
wilcox_observed <- pairwise.wilcox.test(
rich_df$Observed,
rich_df$bmi_group,
p.adjust.method = "fdr"
)
wilcox_observed

# use Wilcoxon rank-sum test for Shannon diversity 
wilcox_shannon <- pairwise.wilcox.test(
rich_df$Shannon,
rich_df$bmi_group,
p.adjust.method = "fdr"
)
wilcox_shannon
savehistory("DOE_dietswap.R")
