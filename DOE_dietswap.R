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

metadata <- data.frame(sample_data(dietswap))
class(metadata)  # Should return "data.frame"

# correction line, in order to tell R that the columns timepoint and timepoint.within.group are categorical variables and not numerical
library(dplyr)
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
freq_plot <- plot_frequencies(x=sample_data(dietswap),
Groups="bmi_group", Factor="nationality")
ggsave("plots/frequency_plot.png", freq_plot, width=8, height=6)

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
rich_plot <- plot_richness(core_diet, x="bmi_group", measures= c("Shannon")) +
geom_boxplot() +
theme_classic() +
theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))
ggsave("plots/richness_plot.png", rich_plot, width=8, height=6)

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
# Create a PCoA plot: principal component analysis allows to visualize and explore patterns that explain our data through multivariate projection 
# odinate the data 
set.seed(100)
ord <- ordinate(physeq = diet_log,
method = "MDS",
distance = "bray")

# prepare eigenvalues to adjust axis
evals <- ord$values$Eigenvalues
# plot PCA
plot_ordination(physeq = diet_log,
ordination = ord,
color = "nationality") +
geom_point(size = 2) +
labs(col = "Nationality") +
coord_fixed(sqrt(evals[2]/evals[1]))

## create a Microbial Abundance Plot: stacked barplots of relative species abundances in each sample group/community
# transform data into relative abundances through the transform-function of the microbiome package
diet_relav <- microbiome::transform(core_diet, "compositional")
#check the transformation process
diet_relav@otu_table@.Data[1:3,1:3]
#filtering data based on label/condition (in our case nationality and bmi group)
afr_lean <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "lean")
afr_over <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "overweight")
afr_obese <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "obese")
aam_lean <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "lean")
aam_over <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "overweight")
aam_obese <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "obese")

#Plot
afr_lean <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "lean")
afr_over <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "overweight")
afr_obese <-  subset_samples(diet_relav, nationality == "AFR" & bmi_group == "obese")
aam_lean <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "lean")
aam_over <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "overweight")
aam_obese <-  subset_samples(diet_relav, nationality == "AAM" & bmi_group == "obese")

plot_composition(afr_lean,
                 taxonomic.level="Genus",
                 average_by="timepoint",
                 otu.sort="abundance",
                 x.label="timepoint")+
  labs(x="Time point",
       y="Abundance",
       title="Native African-Lean")+
  theme(axis.test.x=element_text(angle=0,
                                 hjust=0.5
  ))

# Plot for Native African-Overweight
plot_composition(afr_over,
                 taxonomic.level = "Genus",
                 average_by = "timepoint",
                 otu.sort = "abundance",
                 x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "Native African-Overweight") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# Native African-Lean Plot
afr_lean_plot <- plot_composition(afr_lean,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "Native African-Lean") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# Native African-Overweight Plot
afr_over_plot <- plot_composition(afr_over,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "Native African-Overweight") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# Native African-Obese Plot
afr_obese_plot <- plot_composition(afr_obese,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "Native African-Obese") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# African American-Lean Plot
aam_lean_plot <- plot_composition(aam_lean,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "African American-Lean") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# African American-Overweight Plot
aam_over_plot <- plot_composition(aam_over,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "African American-Overweight") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# African American-Obese Plot
aam_obese_plot <- plot_composition(aam_obese,
                       taxonomic.level = "Genus",
                       average_by = "timepoint",
                       otu.sort = "abundance",
                       x.label = "timepoint") + 
  labs(x = "Time point",
       y = "Abundance",
       title = "African American-Obese") + 
  theme(axis.text.x = element_text(angle = 0, 
                                   hjust = 0.5))

# Save individual plots
ggsave("plots/afr_lean.png", afr_lean_plot, width=8, height=6)
ggsave("plots/afr_overweight.png", afr_over_plot, width=8, height=6)
ggsave("plots/afr_obese.png", afr_obese_plot, width=8, height=6)
ggsave("plots/aam_lean.png", aam_lean_plot, width=8, height=6)
ggsave("plots/aam_overweight.png", aam_over_plot, width=8, height=6)
ggsave("plots/aam_obese.png", aam_obese_plot, width=8, height=6)

# Combine all plots into one figure
combined_plot <- ggarrange(afr_lean_plot, afr_over_plot, afr_obese_plot, 
                         aam_lean_plot, aam_over_plot, aam_obese_plot,
                         ncol=2, nrow=3,
                         common.legend = TRUE)

# Save the combined plot
ggsave("plots/combined_composition_plots.png", 
       combined_plot, 
       width=16, 
       height=18)

