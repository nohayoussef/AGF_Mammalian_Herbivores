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

## creating ordination plots
## PCoA
PCoA_Bray <-ordinate(Physeq, method="PCoA", distance="bray")
P <-plot_ordination(Physeq, PCoA_Bray, "samples", color="Family", shape="Gut")
P=P+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P

PCoA_Unifrac <-ordinate(Physeq, method="PCoA", distance="unifrac")
P2 <-plot_ordination(Physeq, PCoA_Unifrac, "samples", color="Family", shape="Gut")
P2=P2+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P2

PCoA_UnifracW <-ordinate(Physeq, method="PCoA", distance="unifrac", weighted=TRUE)
P3 <-plot_ordination(Physeq, PCoA_UnifracW, "samples", color="Family", shape="Gut")
P3=P3+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P3

## NMDS
NMDS_Bray <-ordinate(Physeq, method="NMDS", distance="bray")
P4 <-plot_ordination(Physeq, NMDS_Bray, "samples", color="Family", shape="Gut")
P4=P4+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P4

NMDS_Unifrac <-ordinate(Physeq, method="NMDS", distance="unifrac")
P5 <-plot_ordination(Physeq, NMDS_Unifrac, "samples", color="Family", shape="Gut")
P5=P5+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P5

NMDS_UnifracW <-ordinate(Physeq, method="NMDS", distance="unifrac", weighted=TRUE)
P6 <-plot_ordination(Physeq, NMDS_UnifracW, "samples", color="Family", shape="Gut")
P6=P6+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P6

##RDA
RDA <-ordinate(Physeq, method="RDA")
P7 <-plot_ordination(Physeq, RDA, "samples", color="Family", shape="Gut")
P7=P7+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
P7

## Distances from Centroids
library(vegan)

## read in the shared file (DM). Column 1 has the sample name and the next columns have the Genera (absolute numbers)
DM <-read.table("DM.txt", header=TRUE, row.names=1)
## Calculate distance matrices based on the community composition
## Bray
Bray <-vegdist(Biom, method="Bray")

## Calculate group centroids for each host factor
## Animal species
Anim_Bray_BD = betadisper(Bray, sample_data(Physeq)$Animal, type="centroid")

## Animal Family
Fam_Bray_BD = betadisper(Bray, sample_data(Physeq)$Family, type="centroid")

## Animal Gut
Gut_Bray_BD = betadisper(Bray, sample_data(Physeq)$Gut, type="centroid")

## Plotting distances to centroids
boxplot(Anim_Bray_BD, ylab="Distance to Centroids")
boxplot(Fam_Bray_BD, ylab="Distance to Centroids")
boxplot(Gut_Bray_BD, ylab="Distance to Centroids")

## creating beta diversity matrices then ordering them
## Factors file contains the variables (animal species, family, gut type, lifestyle)
## Biom is the abundance file but flipped with genera as rows and samples as columns
library(rbiom)
Factors <-read.table("Factors.txt", header=TRUE)
Biom <-read.table("biom.txt", header=TRUE, row.names=1)
Biom <-as.matrix(Biom)

dist_mtx_order = function(d, x){
  m = d %>% as.matrix
  d = as.dist(m[x,x])
  return(d)
}

## reading in host tree, converting it to dist
Host_Sp <-ape::read.tree(file="Host.nwk")
host_tree_d = Host_Sp %>% cophenetic %>% as.dist %>% rescale_dist_mtx

X = labels(host_tree_d)


Bray_o <-dist_mtx_order(Bray, X)
Unifrac_uw <-rbiom::unifrac(Biom, weighted=FALSE, tree=Tree)
Uni_uw_o <-dist_mtx_order(Unifrac_uw, X)
Unifrac_w <-rbiom::unifrac(Biom, weighted=TRUE, tree=Tree)
Uni_w_o <-dist_mtx_order(Unifrac_w, X)

## adonis (PERMANOVA)
Adonis1 <-adonis(Bray_o ~Animal, Factors)
Adonis2 <-adonis(Uni_uw_o ~Animal, Factors)
Adonis3 <-adonis(Uni_w_o ~Animal, Factors)

Adonis4 <-adonis(Bray_o ~Family, Factors)
Adonis5 <-adonis(Uni_uw_o ~Family, Factors)
Adonis6 <-adonis(Uni_w_o ~Family, Factors)

Adonis7 <-adonis(Bray_o ~Lifestyle, Factors)
Adonis8 <-adonis(Uni_uw_o ~Lifestyle, Factors)
Adonis9 <-adonis(Uni_w_o ~Lifestyle, Factors)

Adonis10 <-adonis(Bray_o ~Gut, Factors)
Adonis11 <-adonis(Uni_uw_o ~Gut, Factors)
Adonis12 <-adonis(Uni_w_o ~Gut, Factors)

## The ANOVA table (with R2 and F-test p-value) can be obtained by running the following
Adonis1$aov
Adonis2$aov
Adonis3$aov

Adonis4$aov
Adonis5$aov
Adonis6$aov

Adonis7$aov
Adonis8$aov
Adonis9$aov

Adonis10$aov
Adonis11$aov
Adonis12$aov

## Multivariate methods: MRM, Mantel, Procrustes
Jac <-vegdist(Biom, method="Jaccard")
Jac_o <-dist_mtx_order(Jac, X)

## creating factors distance matrices by Gower transformation
Lifestyle_d <-vegan::vegdist(Factors$Lifestyle, metric="gower")
Lifestyle_d_o <-dist_mtx_order(Lifestyle_d, X)
Gut_d <-vegan::vegdist(Factors$Gut, metric="gower")
Gut_d_o <-dist_mtx_order(Gut_d, X)
Family_d <-vegan::vegdist(Factors$Family, metric="gower")
Family_d_o <-dist_mtx_order(Family_d, X)

Lifestyle_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Gut_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Family_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)

Uni_uw_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Uni_w_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Bray_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Jac_o %>% lapply(function(x) x %>% as.matrix %>% dim)

## Running MRM without subsampling (the second command will display the R2 and p-value)
MRM_Jac <-ecodist::MRM(Jac_o ~ host_tree_d+Family_d_o+Gut_d_o+Lifestyle_d_o)
MRM_Jac
MRM_Bray <-ecodist::MRM(Bray_o ~ host_tree_d+Family_d_o+Gut_d_o+Lifestyle_d_o)
MRM_Bray
MRM_Uni_uw <-ecodist::MRM(Uni_uw_o ~ host_tree_d+Family_d_o+Gut_d_o+Lifestyle_d_o)
MRM_Uni_uw
MRM_Uni_w <-ecodist::MRM(Uni_w_o ~ host_tree_d+Family_d_o+Gut_d_o+Lifestyle_d_o)
MRM_Uni_w

## Running Mantel (the second command will display the coefficient and p-value)
Mantel_Jac <-ecodist::mantel(Jac_o ~ host_tree_d)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Family_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Gut_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Lifestyle_d_o)
Mantel_Jac
Mantel_Bray <-ecodist::mantel(Bray_o ~ host_tree_d)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Family_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Gut_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Lifestyle_d_o)
Mantel_Bray
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ host_tree_d)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Family_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Gut_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Lifestyle_d_o)
Mantel_Uni_uw
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ host_tree_d)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Family_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Gut_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Lifestyle_d_o)
Mantel_Uni_w

## Running Procrustes (the second command will display the sum of squares, correlation coefficient, and p-value)
protest1 <-protest(Jac_o, Gut_d_o, scores="sites")
protest1
protest2 <-protest(Jac_o, Family_d_o, scores="sites")
protest2
protest3 <-protest(Jac_o, host_tree_d, scores="sites")
protest3
protest4 <-protest(Jac_o, Lifestyle_d_o, scores="sites")
protest4
protest5 <-protest(Uni_uw_o, Family_d_o, scores="sites")
protest5
protest6 <-protest(Uni_uw_o, Gut_d_o, scores="sites")
protest6
protest7 <-protest(Uni_uw_o, host_tree_d, scores="sites")
protest7
protest8 <-protest(Uni_uw_o, Lifestyle_d_o, scores="sites")
protest8
protest9 <-protest(Uni_w_o, Gut_d_o, scores="sites")
protest9
protest10 <-protest(Uni_w_o, Family_d_o, scores="sites")
protest10
protest11 <-protest(Uni_w_o, host_tree_d, scores="sites")
protest11
protest12 <-protest(Uni_w_o, Lifestyle_d_o, scores="sites")
protest12
protest13 <-protest(Bray_o, Gut_d_o, scores="sites")
protest13
protest14 <-protest(Bray_o, Family_d_o, scores="sites")
protest14
protest15 <-protest(Bray_o, host_tree_d, scores="sites")
protest15
protest16 <-protest(Bray_o, Lifestyle_d_o, scores="sites")
protest16

## Running MRM with subsampling one per animal species and repeating the analysis 100 times
library(HTSSIP)
# metadata
metadata = Physeq %>%
  phyloseq2df(sample_data) %>%
  dplyr::select(Sample, Animal) %>%
  as.data.frame

rownames(metadata) = metadata$SampleID
metadata[sample_names(Physeq),]
metadata

# number of subsampled trees
ntrees = 100

## randomly selecting one per group
tree_subsample = function(L, df, tree){
  # get subsample (note: subsampling within each species)
  to_keep = df %>% 
    group_by(Animal) %>% 
    sample_n(1) %>%
    .$sample
  # subsampling tree
  to_rm = setdiff(tree$tip.label, to_keep)
  tree = drop.tip(tree, to_rm)
  return(tree)
}

# subsampling trees
df = metadata %>%
  mutate(sample = metadata$Sample) %>%
  dplyr::select(sample, Animal) 
## this following will create a file with 100 trees each with just one from each type of animal
doParallel::registerDoParallel(threads)
host_tree_l = plyr::llply(as.list(1:ntrees), 
                          function(x) tree_subsample(x, df, Host_Sp),
                          .parallel=TRUE)

# tree lengths (should be 100)
# the output of the second line below should be equal to the number of different animal types you have
host_tree_l %>% length %>% print
lapply(host_tree_l, function(x) x$tip.label %>% length) %>% unlist %>% summary

# number of permutated SpecD datasetst
nperm_datasets = 100
# number of permutations per MRM analysis
nperm = 1000
#' randomly selecting one per group
#' L : list of distance matrixes used for MRM
#' df_grps : data.frame (sample, group)
one_per_group = function(L, df_grps, ...){
  # get subsample
  colnames(df_grps) = c('sample', 'Animal')
  df_grps = df_grps %>%
    group_by(Animal) %>%
    sample_n(1)
  # subsetting all matrices
  lapply(L, function(x) dist_mtx_order(x, df_grps$sample))
}

#' MRM on one subsample rep
#' i : rep number
#' L : list of list of distance matrices generated by `one_per_group()`
# nperm : nperm function for MRM
# f : MRM fomulat
mrm_each = function(i, L, f, nperm=99){
  m = L[[i]]
  f = as.formula(f)
  x = ecodist::MRM(f, nperm=nperm, mrank=TRUE)
  # coefficients
  df = x$coef %>% as.data.frame
  colnames(df) = c('coef', 'coef_pval')
  df$variable = rownames(df)
  df$R2 = x$r.squared[1]
  df$pval = x$r.squared[2]
  df$F = x$F.test[1]
  df$rep = i
  return(df)
}

# creating subsample permutations of the distance matrices (starting with Jaccard)
L = list(beta = Jac_o,
         host_phy = host_tree_d,
         Gut = Gut_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Gut+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary

mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()


# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()
plot(p)
write.table(mrm_res, file="mrm_res_Jac")

## below is the same set of commands but repeated by beta diversity measure
# Bray
L = list(beta = Bray_o,
         host_phy = host_tree_d,
         Country = Gut_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Gut+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary

mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()
plot(p)
write.table(mrm_res, file="mrm_res_Bray")

# Unifrac Unweighted
L = list(beta = Uni_uw_o,
         host_phy = host_tree_d,
         Country = Gut_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Gut+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary
mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()
plot(p)
write.table(mrm_res, file="mrm_res_UniUW")

# Weighted Unifrac
L = list(beta = Uni_w_o,
         host_phy = host_tree_d,
         Country = Gut_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Gut+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary

mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           
mrm_res_s

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()
plot(p)
write.table(mrm_res, file="mrm_res_UniW")


 