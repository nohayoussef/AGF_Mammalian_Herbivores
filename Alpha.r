library(readxl)
library(phyloseq)
library(ape)
library(plyr)
library(vegan)
library(rbiom)
library(ggplot2)
options(max.print=999999)

## create phyloseq object with the genera tree. Physeq is an excel file with three sheets:
## Sheet 1: is named OTU and it has the data from the shared file with column 1 header being #OTU ID 
## and containing the genera, Next columns are the samples
## Sheet 2: is named taxon and it has the classification of all genera. Column 1 header being #OTU ID
## and next columns are: Kingdom, Phylum, Class, Order, Genus
## Sheet 3: is named Sample and it has the metadata. Here the first column is named Sample and it has the 
## sample names as rows. Next columns are the metadata entries. Here they are Animal, Family, and GutType.
otu_mat <-read_excel("Physeq.xlsx", sheet="OTU")
tax_mat<- read_excel("Physeq.xlsx", sheet="taxon")
Meta <-read_excel("Physeq.xlsx", sheet="Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
Meta <- Meta  %>%
  tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
Tree <-ape::read.tree(file="Genera_mafft.nwk")
Physeq <-phyloseq(OTU, TAX, samples, Tree)
Physeq

## calculate alpha diversity measures and plot them
Alpha <-estimate_richness(Physeq, measures=c("Observed", "Shannon", "Simpson", "InvSimpson"))
write.table(Alpha, file="Alpha.txt")
p1= plot_richness(Physeq, x="Animal", measures=c("Observed", "Shannon", "Simpson", "InvSimpson"))
p1 + geom_boxplot(data =p1$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
  
p2= plot_richness(Physeq, x="Family", measures=c("Observed", "Shannon", "Simpson", "InvSimpson"))
p2 + geom_boxplot(data =p2$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

p3= plot_richness(Physeq, x="GutType", measures=c("Observed", "Shannon", "Simpson", "InvSimpson"))
p3 + geom_boxplot(data =p3$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
  

## running ANOVA on the alpha diversity results grouped by host factors
## A high F-value and low P-value mean that the variance between groups is significant. 
## The TukeyHSD gives the pairwise p-values for all possible pairs of groups

### on Observed-Animal species
Animal_anova.obs = aov(Alpha$Observed ~ sample_data(Physeq)$Animal)
summary(Animal_anova.obs)
Animal_anova.obs
TukeyHSD(Animal_anova.obs)
### on Observed-Animal Family
Family_anova.obs = aov(Alpha$Observed ~ sample_data(Physeq)$Family)
summary(Family_anova.obs)
Family_anova.obs
TukeyHSD(Family_anova.obs)
### on Observed-GutType
Gut_anova.obs = aov(Alpha$Observed ~ sample_data(Physeq)$GutType)
summary(Gut_anova.obs)
Gut_anova.obs
TukeyHSD(Gut_anova.obs)

### on Shannon-Animal species
Animal_anova.shannon = aov(Alpha$Shannon ~ sample_data(Physeq)$Animal)
summary(Animal_anova.shannon)
Animal_anova.shannon
TukeyHSD(Animal_anova.shannon)
### on Shannon-Animal Family
Family_anova.shannon = aov(Alpha$Shannon ~ sample_data(Physeq)$Family)
summary(Family_anova.shannon)
Family_anova.shannon
TukeyHSD(Family_anova.shannon)
### on Shannon-GutType
Gut_anova.shannon = aov(Alpha$Shannon ~ sample_data(Physeq)$GutType)
summary(Gut_anova.shannon)
Gut_anova.shannon
TukeyHSD(Gut_anova.shannon)

### on Simpson-Animal species
Animal_anova.Simpson = aov(Alpha$Simpson ~ sample_data(Physeq)$Animal)
summary(Animal_anova.Simpson)
Animal_anova.Simpson
TukeyHSD(Animal_anova.Simpson)
### on Simpson-Animal Family
Family_anova.Simpson = aov(Alpha$Simpson ~ sample_data(Physeq)$Family)
summary(Family_anova.Simpson)
Family_anova.Simpson
TukeyHSD(Family_anova.Simpson)
### on Simpson-GutType
Gut_anova.Simpson = aov(Alpha$Simpson ~ sample_data(Physeq)$GutType)
summary(Gut_anova.Simpson)
Gut_anova.Simpson
TukeyHSD(Gut_anova.Simpson)

### on InvSimpson-Animal species
Animal_anova.InvSimpson = aov(Alpha$InvSimpson ~ sample_data(Physeq)$Animal)
summary(Animal_anova.InvSimpson)
Animal_anova.InvSimpson
TukeyHSD(Animal_anova.InvSimpson)
### on InvSimpson-Animal Family
Family_anova.InvSimpson = aov(Alpha$InvSimpson ~ sample_data(Physeq)$Family)
summary(Family_anova.InvSimpson)
Family_anova.InvSimpson
TukeyHSD(Family_anova.InvSimpson)
### on InvSimpson-GutType
Gut_anova.InvSimpson = aov(Alpha$InvSimpson ~ sample_data(Physeq)$GutType)
summary(Gut_anova.InvSimpson)
Gut_anova.InvSimpson
TukeyHSD(Gut_anova.InvSimpson)