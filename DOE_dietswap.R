library(phyloseq)
data("GlobalPatterns")
# View the phyloseq object summary
GlobalPatterns
# View sample metadata (sample characteristics)
sample_data(GlobalPatterns)
# View the OTU table
otu_table(GlobalPatterns)
# View the taxonomy table
tax_table(GlobalPatterns)
# View the phylogenetic tree
phy_tree(GlobalPatterns)
library(phyloseq)
data("DietSwap")
DietSwap
data("dietswap")
data("dietswap")
library(phyloseq)
data("dietswap")
library(microbiome)
data("dietswap")
dietswap
sample_data(dietswap)
data("esophagus")
esophagus
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("microbiomeDataSets")
library(HMP16SData)
BiocManager::install("HMP16SData")
HMP16SData
library(HMP16SData)
HMP16SData
v35_data <- V35()
v35_data
data(bacteria_genus)
library(microbiome)
data(bacteria_genus)
data(bacteria_genus50)
data(enterotype)
enterotype
sample_data(enterotype)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
#upload data dietswaps from package phyloseq (data from mucosa sample analysis regarding different dietary treatments)
data("dietswap")
dietswap
dietswap@otu_table@.Data[1:10,1:7]
sample_data(dietswap)
as.data.frame(sample_data(dietswap))
meta(dietswap)
head(meta(dietswap))
metadata <- metadata %>%
mutate(timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group))
class(metadata)
metadata <- as.data.frame(sample_data(dietswap))
metadata <- metadata %>%
mutate(timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group))
metadata <- metadata %>%
mutate(timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group))
metadata <- as.data.frame(sample_data(dietswap))
library(dplyr)
metadata <- metadata %>%
mutate(timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group))
sample_data(dietswap) <- sample_data(metadata)
metadata <- metadata %>%
mutate(timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group))
metadata <- as.data.frame(sample_data(dietswap))
sample_data(dietswap)
metadata <- data.frame(sample_data(dietswap))
class(metadata)
metadata <- metadata %>%
mutate(
timepoint = as.factor(timepoint),
timepoint.within.group = as.factor(timepoint.within.group)
)
summary(metadata)
diet <- psmelt(dietswap))
diet <- psmelt(dietswap)
plot_frequencies(x=sample_data(dietswap),
Groups="bmi_group", Factor="nationality")
core_diet <- core(x=dietswap,
detection=50,
prevalence=50/100)
diet_log <- transform_sample_counts(core_diet, function(x) log(1+x))
rich <- estimate_richness(core_diet)
rich
rich <- estimate_richness(core_diet, measures = c("Observed", "Shannon"))
rich
plot_richness(core_diet, x="bmi_group", measures= c("Shannon")) +
geom_boxplot() +
theme_classic() +
theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
rich_df <- as.data.frame(rich)
rich_df$bmi_group <- sample_data(core_diet)$bmi_group
wilcox_observed <- pairwise.wilcox.test(
rich_df$Observed,
rich_df$bmi_group,
p.adjust.method = "fdr"
)
wilcox_observed
wilcox_shannon <- pairwise.wilcox.test(
rich_df$Shannon,
rich_df$bmi_group,
p.adjust.method = "fdr"
)
wilcox_shannon
savehistory("DOE_dietswap.R")
