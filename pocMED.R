
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Pocilloporidae samples extracted from all MED run
# 27.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('grid')
library('RColorBrewer')
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

# import normal percent matrix

allShared = read.table("all.7801.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# Import normal taxonomy file from mothur

allTax = read.table('all.7801.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allTax) = allTax[,1]
allTax = allTax[,3:9]
allTax = as.matrix(allTax)

# import meta data

metaFile = read.table('metaDataChem.txt', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:22]

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)

# take a look at the pocilloporid samples

pVerr <- subset_samples(allPhylo, species=='Pocillopora verrucosa')
pDami <- subset_samples(allPhylo, species=='Pocillopora damicornis')
poc <- merge_phyloseq(pVerr, pDami)

sample_data(poc)$names <- factor(sample_names(poc), levels=unique(sample_names(poc)))
sample_data(pVerr)$names <- factor(sample_names(pVerr), levels=unique(sample_names(pVerr)))
sample_data(pDami)$names <- factor(sample_names(pDami), levels=unique(sample_names(pDami)))

pocFilt = filter_taxa(poc, function(x) mean(x) > 0.5, TRUE)

pVerrFilt = filter_taxa(pVerr, function(x) mean(x) > 0.35, TRUE)

# phylum / class bars
# transforming to normal matix and as.factor keeps taxa stacked consistently

taxLevel <- "Genus"

pVerrFiltGlom <- tax_glom(pVerrFilt, taxrank=taxLevel)
pVerrdf <- psmelt(pVerrFiltGlom)

# get total abundance so can make a 'other' column
# had to add ^ and $ characters to make sure grep matches whole word

pVerrdfOther <- pVerrdf

for (j in unique(pVerrdf$Sample)) {
  jFirst = paste('^', j, sep='')
  jBoth = paste(jFirst, '$', sep='')
  rowNumbers = grep(jBoth, pVerrdf$Sample)
  otherValue = 100 - sum(pVerrdf[rowNumbers,"Abundance"])
  newRow = (pVerrdf[rowNumbers,])[1,]
  newRow[,taxLevel] = "other"
  newRow[,"Abundance"] = otherValue
  pVerrdfOther <- rbind(pVerrdfOther, newRow)
}

# need to create my own ggplot colors then replace the last one with gray
# this will ensure that the 'other' category is gray
# need to manually change the tax level after ggplot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggCols <- gg_color_hue(length(unique(pVerrdfOther[,taxLevel])))
ggCols <- head(ggCols, n=-1)

pVerrdfOther$names <- factor(pVerrdfOther$Sample, levels=rownames(metaFile), ordered = TRUE)

theme_set(theme_bw())
ggplot(pVerrdfOther, aes(x=names, y=Abundance, fill=Genus, order = as.factor(Genus))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c(ggCols, "gray")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# SAVE PLOT: EPS 1500 x 600


# let's check what's happening with different Endozoicomonas OTUs (Pocillopora combined)

pocEndo = subset_taxa(poc, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs 

tax_table(pocEndo) <- cbind(tax_table(pocEndo), Strain=taxa_names(pocEndo))
myranks = c("Genus", "Strain")
mylabels = apply(tax_table(pocEndo)[, myranks], 1, paste, sep="", collapse="_")
tax_table(pocEndo) <- cbind(tax_table(pocEndo), catglab=mylabels)

sample_data(pocEndo)$names <- factor(sample_names(pocEndo), levels=unique(sample_names(pocEndo)))

pocEndoFilt = filter_taxa(pocEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(pocEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_manual(values=cols2) +
  facet_grid(~site, scales='free', space='free_x')

# let's check what's happening with different Endozoicomonas OTUs (Pocillopora verrucosa)

pVerrEndo = subset_taxa(pVerr, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs 

tax_table(pVerrEndo) <- cbind(tax_table(pVerrEndo), Strain=taxa_names(pVerrEndo))
myranks = c("Genus", "Strain")
mylabels = apply(tax_table(pVerrEndo)[, myranks], 1, paste, sep="", collapse="_")
tax_table(pVerrEndo) <- cbind(tax_table(pVerrEndo), catglab=mylabels)

sample_data(pVerrEndo)$names <- factor(sample_names(pVerrEndo), levels=unique(sample_names(pVerrEndo)))

pVerrEndoFilt = filter_taxa(pVerrEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(pVerrEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Set1') +
  facet_grid(~site, scales='free', space='free_x')

# SAVE EPS 1500 x 700

# ordination for the pocillopora samples

theme_set(theme_bw())
pocOrd <- ordinate(poc, "NMDS", "bray")
pocOrdPlot <- plot_ordination(poc, pocOrd, type='samples', color='site', title='poc') +
  scale_color_manual(values=cols) 

# not many patterns apparent
# subset based on poc type and re-draw

pocType3 <- subset_samples(poc, pocType=='type3')
pocType3Ord <- ordinate(pocType3, "NMDS", "bray")
plot_ordination(pocType3, pocType3Ord, type = 'samples', color='site', title='pocType3') +
  scale_color_manual(values=cols) 

pocTypeUnknown <- subset_samples(poc, pocType=='type?')
pocTypeUnknownOrd <- ordinate(pocTypeUnknown, "NMDS", "bray")
plot_ordination(pocTypeUnknown, pocTypeUnknownOrd, type = 'samples', color='site', title='pocTypeUnknown') +
  scale_color_manual(values=cols) 

pocType5 <- subset_samples(poc, pocType=='type5')
pocType5Ord <- ordinate(pocType5, "NMDS", "bray")
plot_ordination(pocType5, pocType5Ord, type = 'samples', color='site', title='pocType5') +
  scale_color_manual(values=cols) 

pocType1 <- subset_samples(poc, pocType=='type1')

pocType5andUnknown <- merge_phyloseq(pocTypeUnknown, pocType5, pocType1)
pocType5andUnknownOrd <- ordinate(pocType5andUnknown, "NMDS", "bray")
plot_ordination(pocType5andUnknown, pocType5andUnknownOrd, type = 'samples', color='site', title='pocType5andUnknown', label='names') +
  scale_color_manual(values=cols) 

# SAVE AS 700 x 532

# check what happens with Endos across poc Types

pocEndoFilt = filter_taxa(pocEndo, function(x) mean(x) > 0.1, TRUE)

theme_set(theme_bw())
plot_bar(pocEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  #scale_fill_brewer(type='qual', palette = 'Set1') +
  facet_grid(~pocType, scales='free', space='free_x')


# ordination for the pocillopora verrucosa

theme_set(theme_bw())
pVerrOrd <- ordinate(pVerr, "NMDS", "bray")
plot_ordination(pVerr, pVerrOrd, type = 'samples', color='site', title='pVerr')

# overlay the best Spearman correlations from mothur

pVerrMothurSpearman <- c("MED000007996","MED000007267","MED000007994","MED000007995","MED000007271","MED000007997","MED000007227","MED000006335","MED000003910","MED000000074")


arrowmatrix <- data.frame(labels = pVerrMothurSpearman)
rownames(arrowmatrix) <- arrowmatrix[,1]
arrowdf <- data.frame(arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Order"], arrowmatrix)
arrowdf <- cbind(arrowdf, pVerrOrd$species[pVerrMothurSpearman, ])

# also combine MED node with genus taxonomy for labelling

mylabels = apply(arrowdf[, c("labels", "labels.1")], 1, paste, sep="", collapse="_")
arrowdf <- cbind(arrowdf, catglab=mylabels)

# now create labels, arrows, etc.

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = catglab, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))

pVerrFiltPlot <- plot_ordination(pVerr, pVerrOrd, type = 'samples', color='site', title='pVerr')

pVerrFiltPlotArrow <- pVerrFiltPlot + 
  geom_segment(arrowmap, size = 1, data = arrowdf, color = "grey",  arrow = arrowhead, show_guide = FALSE) +
  scale_color_manual(values=cols) +
  geom_text(labelmap, size = 2, data = arrowdf)

pVerrFiltPlotArrow



####################################################################
### fit chemistry to the coral microbiomes
####################################################################

waterQual <- c("temp", "salinity", "Domg", "pH")
nutrients <- c("PO4", "N.N", "silicate", "NO2", "NH4")
FCM <- c("prok", "syn", "peuk", "pe.peuk", "Hbact")

chemNoNA <- na.omit(metaFile[sample_names(poc),FCM])
pocNoNA <- prune_samples(rownames(chemNoNA), poc)

sample_names(poc)
sample_names(pocNoNA)

theme_set(theme_bw())
pocOrdNoNA <- ordinate(pocNoNA, "NMDS", "bray")
pocOrdNoNAPlot <- plot_ordination(pocNoNA, pocOrdNoNA, type = 'samples', color='site', title='poc') +
  geom_point(size=3) +
  scale_color_manual(values=c(cols)) 
pocOrdNoNAPlot

pointsNoNA <- pocOrdNoNA$points[rownames(chemNoNA),]

chemFit <- envfit(pointsNoNA, env = chemNoNA, na.rm=TRUE)

chemFit.scores <- as.data.frame(scores(chemFit, display= "vectors"))
chemFit.scores <- cbind(chemFit.scores, Species = rownames(chemFit.scores))

# create arrow info again

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = rownames(chemFit.scores))
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, size=1.5, label = rownames(chemFit.scores))
arrowhead = arrow(length = unit(0.25, "cm"))

pocOrdNoNAPlot + 
  coord_fixed() +
  geom_segment(arrowmap, size = 0.5, data = chemFit.scores, color = "black",  arrow = arrowhead, show_guide = FALSE) +
  geom_text(labelmap, size = 3, data = chemFit.scores)

# SAVE AS 700 x 532

***VECTORS

MDS1     MDS2     r2 Pr(>r)   
temp     -0.98044  0.19681 0.4273  0.002 **
  salinity  0.84477 -0.53513 0.3543  0.002 **
  Domg      0.27923 -0.96022 0.0077  0.911   
pH        0.62615 -0.77971 0.2928  0.014 * 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.

MDS1     MDS2     r2 Pr(>r)  
PO4      -0.05164  0.99867 0.0542  0.411  
N.N      -0.02462  0.99970 0.0612  0.368  
silicate -0.18433 -0.98286 0.0435  0.493  
NO2       0.41953 -0.90774 0.1097  0.157  
NH4      -0.50160 -0.86510 0.1784  0.036 *

  MDS1     MDS2     r2 Pr(>r)  
prok    -0.55526  0.83168 0.2551  0.011 *
  syn     -0.88714  0.46151 0.1525  0.052 .
peuk    -0.92447  0.38126 0.1080  0.150  
pe.peuk -0.90889  0.41703 0.0684  0.332  
Hbact   -0.79815  0.60247 0.1914  0.026 *
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.


####################################################################
### SIMPROF analysis to check which samples fall into 'groups' without any *a priori* assumptions
####################################################################

# Need to import the shared file containing just poc OTUs, then calcualte the simprof clusters based on the braycurtis metric. 

pVerrShared = otu_table(pVerr)
class(pVerrShared) <- "numeric"

pVerrSIMPROF <- simprof(pVerrShared, num.expected=1000, num.simulated=999, method.cluster='average', method.distance='braycurtis', method.transform='squareroot', alpha=0.05, sample.orientation='row', silent=FALSE)

simprof.plot(pVerrSIMPROF, leafcolors=NA, plot=TRUE, fill=TRUE, leaflab="perpendicular", siglinetype=1)

# SAVE EPS 1500 x 700

# I'll try and overlay the significant clusters on top of the nMDS. 
# After calculating the clusters, make a data frame of the results and add to previous nMDS plot. Need to add these groups to the nMDS data.frame - I'll do a loop for this.

simprofCLUSTERS = data.frame()

for(j in 1:length(pocSIMPROF$significantclusters)){
  if(length(pocSIMPROF$significantclusters[[j]]) > 2){
    simprofCLUSTERS <- rbind(simprofCLUSTERS, cbind(j, pocSIMPROF$significantclusters[[j]]))
  }
}

rownames(simprofCLUSTERS) <- simprofCLUSTERS[,2]
colnames(simprofCLUSTERS) <- c("simprofCLUSTER", "group")
pocFiltDF <- as.data.frame(pocOrdPlot$data)
metaSIMP <- merge(pocFiltDF, simprofCLUSTERS, by="row.names")

#plot over nMDS

df_ell <- data.frame()
for(g in levels(metaSIMP$simprofCLUSTER)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metaSIMP[metaSIMP$simprofCLUSTER==g,], ellipse(cor(NMDS1, NMDS2), scale=c(sd(NMDS1), sd(NMDS2)), centre=c(mean(NMDS1), mean(NMDS2))))),site=g))
}

nMDSsimprof <- ggplot(data=metaSIMP, aes(x=NMDS1, y=NMDS2, color=site)) +
  geom_point() +
  scale_color_hue(limits=levels(droplevels(metaSIMP$site)))
nMDSsimprof

nMDS.mean <- aggregate(df_ell, list(df_ell$site), FUN=mean)

#write.table(df_ell, file="simpBrayIdent.txt")
#write.table(nMDS.mean, file="simpBrayIdentMean.txt")

#df_ell_saved <- (read.table("simpBrayIdent.txt"))
#df_ell_saved$site <- factor(df_ell_saved$site)
#nMDS.mean.saved <- read.table("simpBrayIdentMean.txt")

#unique(df_ell$site)

#wanted <- c("2", "3",  "5",  "8",  "10", "16", "21", "26", "27")

#df_ell2 <- subset(df_ell, site==wanted)

pocOrdPlot +
  geom_path(data=df_ell, aes(x=x, y=y), color='black', size=0.5, linetype=2, show_guide=FALSE) +
  annotate("text", x=nMDS.mean$x, y=nMDS.mean$y, label=nMDS.mean$Group.1)

# SAVE AS 700 x 532


####################################################################
### do adonis / PERMANOVA sig testing 1.4.15
####################################################################

# use an NMDS distance matrix for significance testing

pVerrDist <-  phyloseq::distance(pVerr, "bray")
pVerrNMDS <- ordinate(pVerr, "NMDS", pVerrDist)

pVerrADONIS <- adonis(pVerrDist ~ site, as(sample_data(pVerr), "data.frame"))

adonis(formula = pVerrDist ~ site, data = as(sample_data(pVerr),      "data.frame")) 

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
site       3    2.6398 0.87995    2.31 0.13093  0.001 ***
  Residuals 46   17.5225 0.38092         0.86907           
Total     49   20.1624                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# could also check pairwise using a reduced dataset
# Red Sea vs Indonesia

pVerrRed <- subset_samples(pVerr, site=="RedSea")
pVerrIndon <- subset_samples(pVerr, site=="Indonesia")
pVerrRedIndon <- merge_phyloseq(pVerrRed, pVerrIndon)

pVerrRedIndonDist <-  phyloseq::distance(pVerrRedIndon, "bray")
pVerrRedIndonNMDS <- ordinate(pVerrRedIndon, "NMDS", pVerrRedIndonDist)

adonis(pVerrRedIndonDist ~ site, as(sample_data(pVerrRedIndon), "data.frame"))

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
site       1    1.5204 1.52043  4.1981 0.09948  0.001 ***
  Residuals 38   13.7626 0.36217         0.90052           
Total     39   15.2830                 1.00000   

# Indo vs Micro

pVerrMicro <- subset_samples(pVerr, site=="Micronesia")
pVerrIndon <- subset_samples(pVerr, site=="Indonesia")
pVerrMicroIndon <- merge_phyloseq(pVerrMicro, pVerrIndon)

pVerrMicroIndonDist <-  phyloseq::distance(pVerrMicroIndon, "bray")
pVerrMicroIndonNMDS <- ordinate(pVerrMicroIndon, "NMDS", pVerrMicroIndonDist)

adonis(pVerrMicroIndonDist ~ site, as(sample_data(pVerrMicroIndon), "data.frame"))

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
site       1     0.738 0.73804  2.0101 0.06089  0.013 *
  Residuals 31    11.382 0.36717         0.93911         
Total     32    12.120                 1.00000 

# Mal vs Micro

pVerrMicro <- subset_samples(pVerr, site=="Micronesia")
pVerrMal <- subset_samples(pVerr, site=="Maldives")
pVerrMicroMal <- merge_phyloseq(pVerrMicro, pVerrMal)

pVerrMicroMalDist <-  phyloseq::distance(pVerrMicroMal, "bray")
pVerrMicroMalNMDS <- ordinate(pVerrMicroMal, "NMDS", pVerrMicroMalDist)

adonis(pVerrMicroMalDist ~ site, as(sample_data(pVerrMicroMal), "data.frame"))

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
site       1    0.3710 0.37103 0.78944 0.08982  0.951
Residuals  8    3.7599 0.46999         0.91018       
Total      9    4.1310                 1.00000 


## ok now to check which are different using Tukey tests

# calculate multivariate dispersions

siteFactor <- as(sample_data(pVerr), "data.frame")

mod <- betadisper(pVerrDist, siteFactor$site)
mod

# perform test

anova(mod)

Analysis of Variance Table

Response: Distances
Df  Sum Sq   Mean Sq F value Pr(>F)
Groups     3 0.01803 0.0060109   0.639 0.5938
Residuals 46 0.43270 0.0094066                     
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# permutation test for F

permutest(mod, pairwise = TRUE)

Pairwise comparisons:
  (Observed p-value below diagonal, permuted p-value above diagonal)
            Indonesia Maldives Micronesia RedSea
Indonesia             0.71900    0.20800  0.468
Maldives     0.71090             0.12200  0.920
Micronesia   0.18456  0.11800             0.556
RedSea       0.47871  0.91891    0.52935        

# Tukey's honest significance differences

mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)

