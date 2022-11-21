library(NST)
## read in the shared file (DM). Column 1 has the sample name and the next columns have the Genera (absolute numbers)

DM <-read.table("DM.txt", header=TRUE, row.names=1)
## read in files with the metadata. I have them broken down into separate files:
## one for Animal species, one for family, one for gut type, and so on. These files are two columns:
## sample ID in column 1 and metadat in column 2
Animal <-read.table("Animal.txt", row.names=1, header=TRUE)
Family <-read.table("Family.txt", row.names=1, header=TRUE)
Gut <-read.table("Gut.txt", row.names=1, header=TRUE)

## For each of these factors, calculate tNST based on Jaccard index
tnstoutAnim_Jac=NST::tNST(comm=DM, dist.method="jaccard", group=Animal,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutAnim_Jac$index.grp, file = "tnstindex.Animal.jaccard.txt")
write.table(tnstoutAnim_Jac$index.pair, file = "tnstindexpair.Animal.jaccard.txt")
write.table(tnstoutAnim_Jac$index.pair.grp, file = "tnstindexpairgrp.Animal.jaccard.txt")

tnstoutFamily_Jac=NST::tNST(comm=DM, dist.method="jaccard", group=Family,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutFamily_Jac$index.grp, file = "tnstindex.Family.jaccard.txt")
write.table(tnstoutFamily_Jac$index.pair, file = "tnstindexpair.Family.jaccard.txt")
write.table(tnstoutFamily_Jac$index.pair.grp, file = "tnstindexpairgrp.Family.jaccard.txt")

tnstoutGut_Jac=NST::tNST(comm=DM, dist.method="jaccard", group=Gut,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutGut_Jac$index.grp, file = "tnstindex.Gut.jaccard.txt")
write.table(tnstoutGut_Jac$index.pair, file = "tnstindexpair.Gut.jaccard.txt")
write.table(tnstoutGut_Jac$index.pair.grp, file = "tnstindexpairgrp.Gut.jaccard.txt")

## For each of these factors, bootstrap tNST based on Jaccard index
tnst.btAnim.Jac=NST::nst.boot(nst.result=tnstoutAnim_Jac, group=Animal, nworker=32, out.detail=TRUE)
write.table(tnst.btAnimal.Jac$summary,file = "tnstboot.Animal.Jac.txt")
write.table(tnst.btAnimal.Jac$detail$NST.boot,file = "tnstboot.Animal.Jac.detail.txt")
tnst.btFamily.Jac=NST::nst.boot(nst.result=tnstoutFamily_Jac, group=Family, nworker=32, out.detail=TRUE)
write.table(tnst.btFamily.Jac$summary,file = "tnstboot.Family.Jac.txt")
write.table(tnst.btFamily.Jac$detail$NST.boot,file = "tnstboot.Family.Jac.detail.txt")
tnst.btGut.Jac=NST::nst.boot(nst.result=tnstoutGut_Jac, group=Gut, nworker=32, out.detail=TRUE)
write.table(tnst.btGut.Jac$summary,file = "tnstboot.Gut.Jac.txt")
write.table(tnst.btGut.Jac$detail$NST.boot,file = "tnstboot.Gut.Jac.detail.txt")

## For each of these factors, calculate tNST based on Bray index
tnstoutAnim_Bray=NST::tNST(comm=DM, dist.method="Bray", group=Animal,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutAnim_Bray$index.grp, file = "tnstindex.Animal.Bray.txt")
write.table(tnstoutAnim_Bray$index.pair, file = "tnstindexpair.Animal.Bray.txt")
write.table(tnstoutAnim_Bray$index.pair.grp, file = "tnstindexpairgrp.Animal.Bray.txt")

tnstoutFamily_Bray=NST::tNST(comm=DM, dist.method="Bray", group=Family,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutFamily_Bray$index.grp, file = "tnstindex.Family.Bray.txt")
write.table(tnstoutFamily_Bray$index.pair, file = "tnstindexpair.Family.Bray.txt")
write.table(tnstoutFamily_Bray$index.pair.grp, file = "tnstindexpairgrp.Family.Bray.txt")

tnstoutGut_Bray=NST::tNST(comm=DM, dist.method="Bray", group=Gut,
                     abundance.weighted=FALSE, nworker=32,
                     null.model="PF", output.rand = TRUE,
                     SES = TRUE, RC = TRUE)
write.table(tnstoutGut_Bray$index.grp, file = "tnstindex.Gut.Bray.txt")
write.table(tnstoutGut_Bray$index.pair, file = "tnstindexpair.Gut.Bray.txt")
write.table(tnstoutGut_Bray$index.pair.grp, file = "tnstindexpairgrp.Gut.Bray.txt")

## For each of these factors, bootstrap tNST based on Bray index
tnst.btAnim.Bray=NST::nst.boot(nst.result=tnstoutAnim_Bray, group=Animal, nworker=32, out.detail=TRUE)
write.table(tnst.btAnimal.Bray$summary,file = "tnstboot.Animal.Bray.txt")
write.table(tnst.btAnimal.Bray$detail$NST.boot,file = "tnstboot.Animal.Bray.detail.txt")
tnst.btFamily.Bray=NST::nst.boot(nst.result=tnstoutFamily_Bray, group=Family, nworker=32, out.detail=TRUE)
write.table(tnst.btFamily.Bray$summary,file = "tnstboot.Family.Bray.txt")
write.table(tnst.btFamily.Bray$detail$NST.boot,file = "tnstboot.Family.Bray.detail.txt")
tnst.btGut.Bray=NST::nst.boot(nst.result=tnstoutGut_Bray, group=Gut, nworker=32, out.detail=TRUE)
write.table(tnst.btGut.Bray$summary,file = "tnstboot.Gut.Bray.txt")
write.table(tnst.btGut.Bray$detail$NST.boot,file = "tnstboot.Gut.Bray.detail.txt")

library(ggplot2)
## plotting bootstarpping results
Animal_tnst_bray_bt <-read.table("tnstboot.Anim.bray.detail.txt", header=TRUE)
Family_tnst_bray_bt <-read.table("tnstboot.Family.bray.detail.txt", header=TRUE)
Gut_tnst_bray_bt <-read.table("tnstboot.Gut.bray.detail.txt", header=TRUE)
Animal_tnst_Jac_bt <-read.table("tnstboot.Anim.Jac.detail.txt", header=TRUE)
Family_tnst_Jac_bt <-read.table("tnstboot.Family.Jac.detail.txt", header=TRUE)
Gut_tnst_Jac_bt <-read.table("tnstboot.Gut.Jac.detail.txt", header=TRUE)
AnimalOrdBray <-factor(Animal_tnst_bray_bt$Animal, levels=c("Cow", "Goat", "Sheep", "Deer", "Horse"))
AnimalOrdJac <-factor(Animal_tnst_Jac_bt$Animal, levels=c("Cow", "Goat", "Sheep", "Deer", "Horse"))
FamilyOrdBray <-factor(Family_tnst_bray_bt$Family, levels=c("Bovidae", "Cervidae", "Camelidae", "Equidae", "Elephantidae"))
FamilyOrdJac <-factor(Family_tnst_Jac_bt$Family, levels=c("Bovidae", "Cervidae", "Camelidae", "Equidae", "Elephantidae"))
GutOrdBray <-factor(Gut_tnst_bray_bt$Gut, levels=c("Foregut", "Pseudoruminant", "Hindgut"))
GutOrdJac <-factor(Gut_tnst_Jac_bt$Gut, levels=c("Foregut", "Pseudoruminant", "Hindgut"))
boxplot(tNST~AnimalOrdBray, data=Animal_tnst_bray_bt, ylab="tNST", las=2, col="grey")
boxplot(tNST~AnimalOrdJac, data=Animal_tnst_Jac_bt, ylab="tNST", las=2, col="grey")
boxplot(tNST~FamilyOrdBray, data=Family_tnst_bray_bt, ylab="tNST", las=2, col="grey")
boxplot(tNST~FamilyOrdJac, data=Family_tnst_Jac_bt, ylab="tNST", las=2, col="grey")
boxplot(tNST~GutOrdBray, data=Gut_tnst_bray_bt, ylab="tNST", las=2, col="grey")
boxplot(tNST~GutOrdJac, data=Gut_tnst_Jac_bt, ylab="tNST", las=2, col="grey")


## calculate values of beta net relatedness index (B-NRI), 
## and modified Raup-Crick metric based on Bray Curtis metric (RCBray) 
library(iCAMP)
library(ape)
library(vegan)
library(reshape2)

Tree <-ape::read.tree(file="Genera_Unrooted.nwk")
Tree_d <-cophenetic.phylo(Tree)
save.wd="/projects/cas002/R"
setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = Tree, wd=save.wd, nworker = 32, memory.G = 50)
}else{
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}
## calculate bNRI and RCbray
AllNRI=bNRIn.p(DM, Tree_d, nworker = 32, memo.size.GB = 50, weighted = TRUE, rand = 1000)
AllNRI_RC=bNRIn.p(DM, Tree_d, nworker = 32, memo.size.GB = 50, weighted = TRUE, rand = 1000, sig.index="RC")
BNRI_2c=melt(AllNRI$index)
write.table(BNRI_2c, file="BNRI_All")
BNRI_RC_2c=melt(AllNRI_RC$index)
write.table(BNRI_RC_2c, file="BNRI_All_RC")

## From the values in the above two files, Excel was used: Values of bNRI were used first to partition selective processes into 
## homogenous (number of pairwise comparisons with b-NRI values < -2), and heterogenous selection 
## (number of pairwise comparisons with b-NRI values > 2). 
## All other pairwise comparisons (with absolute b-NRI values < 2) are considered contributing to stochastic processes 
## (not assigned to selection), and can be further broken down into dispersal and drift based on the taxonomic diversity (values of RCBray). 
## the number of pairwise comparisons with absolute values of RCBray < 0.95 are considered contributing to drift, 
## while the number of pairwise comparisons with absolute values of RCBray > 0.95 are considered contributing to dispersal. 
## This last fraction can be further broken down into homogenizing dispersal (RCBray values <-0.95), and dispersal limitation (RCBray values >0.95). 
## The contribution of each of these processes (homogenous selection, heterogenous selection, homogenizing dispersal, dispersal limitation, and drift) to the total 
## community assembly was calculated from the corresponding number of pairwise comparisons falling into each category as a percentage of all pairwise comparisons.
