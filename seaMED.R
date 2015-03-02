
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# seawater samples extracted from #7 all MED run
# 28.1.2015

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

metaFile = read.table('metaData2.test.txt', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]

### Create phyloseq object

OTU = otu_table(allShared, taxa_are_rows = FALSE)
TAX = tax_table(allTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)

# bar chart of s.pistillata phyla or class - the names parameter preserves order

sea <- subset_samples(allPhylo, species=='seawater')
sample_data(sea)$names <- factor(sample_names(sea), levels=rownames(metaFile), ordered = TRUE)
seaFilt = filter_taxa(sea, function(x) mean(x) > 0.1, TRUE)

taxLevel <- "Class"

seaFiltGlom <- tax_glom(seaFilt, taxrank=taxLevel)
physeqdf <- psmelt(seaFiltGlom)

# get total abundance so can make a 'other' column
# had to add ^ and $ characters to make sure grep matches whole word

physeqdfOther <- physeqdf

for (j in unique(physeqdf$Sample)) {
  jFirst = paste('^', j, sep='')
  jBoth = paste(jFirst, '$', sep='')
  rowNumbers = grep(jBoth, physeqdf$Sample)
  otherValue = 100 - sum(physeqdf[rowNumbers,"Abundance"])
  newRow = (physeqdf[rowNumbers,])[1,]
  newRow[,taxLevel] = "other"
  newRow[,"Abundance"] = otherValue
  physeqdfOther <- rbind(physeqdfOther, newRow)
}

# need to create my own ggplot colors then replace the last one with gray
# this will ensure that the 'other' category is gray

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggCols <- gg_color_hue(length(unique(physeqdfOther[,taxLevel])))
ggCols <- head(ggCols, n=-1)

physeqdfOther$names <- factor(physeqdfOther$Sample, levels=rownames(metaFile), ordered = TRUE)

theme_set(theme_bw())
ggplot(physeqdfOther, aes(x=names, y=Abundance, fill=Class, order = as.factor(Class))) +
  geom_bar(stat="identity", colour="black") +
  scale_fill_manual(values=c(ggCols, "gray")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  facet_grid(~site, scales='free', space='free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# SAVE EPS 1500 x 600

# let's check what's happening with different Endozoicomonas OTUs

seaEndo = subset_taxa(sea, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs

tax_table(seaEndo) <- cbind(tax_table(seaEndo), Strain=taxa_names(seaEndo))

myranks = c("Genus", "Strain")
mylabelssea = apply(tax_table(seaEndo)[, myranks], 1, paste, sep="", collapse="_")
tax_table(seaEndo) <- cbind(tax_table(seaEndo), catglab=mylabelssea)

sample_data(seaEndo)$names <- factor(sample_names(seaEndo), levels=rownames(metaFile), ordered = TRUE)

seaEndoFilt = filter_taxa(seaEndo, function(x) mean(x) > 0.0001, TRUE)

theme_set(theme_bw())
plot_bar(seaEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.6)) +
  #scale_fill_brewer(type='qual', palette = 'Dark2') +
  facet_grid(~site, scales='free', space='free_x')


# ordination for the seawater samples

theme_set(theme_bw())
seaOrd <- ordinate(sea, "NMDS", "bray")
seaOrdPlot <- plot_ordination(sea, seaOrd, type = 'samples', color='site', title='seawater') +
  geom_point(size=3) +
  scale_color_manual(values=c(cols)) 
seaOrdPlot

# let's see if the chemistry is correlated with this
# first do PCA of just chemical data - can do this using R base packages
# extract chemistry data from data frame


nutFile = read.table('Nutrients.txt', header=T, sep='\t')
rownames(nutFile) = nutFile[,1]
nutFile2 = nutFile[,3:(length(nutFile))]

nutFile = na.omit(read.table('FCMeditedSites.txt', header=T, sep='\t'))
rownames(nutFile) = nutFile[,1]
nutFile2 = nutFile[,3:(length(nutFile))]

nutFile = read.table('waterQual.txt', header=T, sep='\t')
nutFile = nutFile[c("sample", "site", "temp", "salinity", "Domg", "pH")]
rownames(nutFile) = nutFile[,1]
nutFile2 = nutFile[,3:(length(nutFile))]

PCA.res <- prcomp(na.omit(nutFile2), scale.=T, center=TRUE)
biplot(PCA.res)
summary(PCA.res)

# better to plot this in ggplot
# extract x (sample location) and rotation (arrow location)

nutPCAggdata <- data.frame(PCA.res$x)
nutPCArotation <- data.frame(PCA.res$rotation)

# get correct multiplier value for arrows

multi <- vegan:::ordiArrowMul(nutPCArotation)
nutPCAarrows <- multi*nutPCArotation

# fix site coloring

#site <- c("MaggieIs", "MaggieIs", "MaggieIs", "Ningaloo", "Ningaloo", "Ningaloo", "Ningaloo", "AmericanSamoa", "AmericanSamoa", "Maldives", "Maldives", "Maldives", "Maldives", "Maldives", "Maldives", "RedSea", "RedSea", "RedSea", "RedSea", "RedSea", "RedSea")
#nutPCAggdata <- cbind(nutPCAggdata, site)

nutPCAggdata <- cbind(nutPCAggdata, nutFile[2])

# create arrow info

arrowmap <- aes(xend = PC1, yend = PC2, x = 0, y = 0, shape = NULL, color = NULL, label = rownames(nutPCAarrows))
labelmap <- aes(x = PC1, y = PC2 + 0.04, shape = NULL, color = NULL, size=1.5, label = rownames(nutPCAarrows))
arrowhead = arrow(length = unit(0.02, "npc"))

nutPCA <- ggplot(nutPCAggdata, title='PCA of the nutrient data') +
  geom_point(aes(x=PC1, y=PC2, color=site), size=3) +
  scale_color_manual(values=c(cols)) 
nutPCA 

nutPCA + geom_segment(arrowmap, size = 0.5, data = nutPCAarrows, color = "black",  arrow = arrowhead, show_guide = FALSE) +
  geom_text(labelmap, size = 3, data = nutPCAarrows)

# SAVE AS 700 x 532

# nutrients
Importance of components:
  PC1    PC2    PC3    PC4    PC5
Standard deviation     1.5436 0.9727 0.9093 0.8045 0.4438
Proportion of Variance 0.4766 0.1893 0.1654 0.1294 0.0394
Cumulative Proportion  0.4766 0.6658 0.8312 0.9606 1.0000

# FCM
Importance of components:
  PC1    PC2    PC3     PC4     PC5
Standard deviation     1.8279 0.9107 0.8127 0.32416 0.25255
Proportion of Variance 0.6682 0.1659 0.1321 0.02102 0.01276
Cumulative Proportion  0.6682 0.8341 0.9662 0.98724 1.00000

# water qual
Importance of components:
  PC1    PC2     PC3     PC4
Standard deviation     1.6603 0.8550 0.59512 0.39790
Proportion of Variance 0.6891 0.1827 0.08854 0.03958
Cumulative Proportion  0.6891 0.8719 0.96042 1.00000


#########################################################
# fit chemistry data to my seawater bacteria ordinations
#########################################################

waterQual <- c("temp", "salinity", "Domg", "pH")
nutrients <- c("PO4", "N.N", "silicate", "N02", "NH4")
FCM <- c("prok", "syn", "peuk", "pe.peuk", "Hbact")

chemNoNA <- na.omit(metaFile[sample_names(sea),nutrients])
seaNoNA <- prune_samples(rownames(chemNoNA), sea)

sample_names(sea)
sample_names(seaNoNA)

theme_set(theme_bw())
seaOrdNoNA <- ordinate(seaNoNA, "NMDS", "bray")
seaOrdNoNAPlot <- plot_ordination(seaNoNA, seaOrdNoNA, type = 'samples', color='site', title='seawater') +
  geom_point(size=3) +
  scale_color_manual(values=c(cols)) 
seaOrdNoNAPlot

pointsNoNA <- seaOrdNoNA$points[rownames(chemNoNA),]

chemFit <- envfit(pointsNoNA, env = chemNoNA, na.rm=TRUE)

chemFit.scores <- as.data.frame(scores(chemFit, display= "vectors"))
chemFit.scores <- cbind(chemFit.scores, Species = rownames(chemFit.scores))

# create arrow info again

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = rownames(chemFit.scores))
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, size=1.5, label = rownames(chemFit.scores))
arrowhead = arrow(length = unit(0.25, "cm"))

seaOrdNoNAPlot + 
  coord_fixed() +
  geom_segment(arrowmap, size = 0.5, data = chemFit.scores, color = "black",  arrow = arrowhead, show_guide = FALSE) +
  geom_text(labelmap, size = 3, data = chemFit.scores)

# SAVE AS 700 x 532

        MDS1      MDS2     r2 Pr(>r)   
prok    -0.164138  0.986440 0.1063  0.099 . 
syn      0.202584 -0.979260 0.2228  0.007 **
peuk     0.273535 -0.961860 0.2637  0.006 **
pe.peuk  0.020191 -0.999800 0.0048  0.907   
Hbact    0.128122 -0.991760 0.0986  0.130   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.

        MDS1     MDS2     r2 Pr(>r)    
temp      0.95414  0.29937 0.5323  0.001 ***
salinity -0.95573 -0.29424 0.5172  0.001 ***
Domg     -0.42780 -0.90387 0.1527  0.092 .  
pH       -0.84785 -0.53023 0.4743  0.002 ** 
  ---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.

MDS1     MDS2     r2 Pr(>r)    
PO4       0.90349 -0.42861 0.4581  0.001 ***
N.N       0.86178 -0.50728 0.2030  0.016 *  
silicate  0.48720 -0.87329 0.3933  0.002 ** 
N02       0.68859 -0.72516 0.7457  0.001 ***
NH4       0.88976 -0.45643 0.1836  0.036 *  
  ---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.




#########################################################
# can also use Mothur to get Pearson correlations
#########################################################

# sample_names(sea)

seaMothurSpearman <- c("MED000007608","MED000006893","MED000005444","MED000004065","MED000007963","MED000006897","MED000007201","MED000004072","MED000006014","MED000006278")

arrowmatrix <- data.frame(labels = seaMothurSpearman)
rownames(arrowmatrix) <- arrowmatrix[,1]
arrowdf <- data.frame(arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Class"], arrowmatrix)
arrowdf <- cbind(arrowdf, seaOrd$species[seaMothurSpearman, ])

# also combine MED node with genus taxonomy for labelling

mylabels = apply(arrowdf[, c("labels", "labels.1")], 1, paste, sep="", collapse="_")
arrowdf <- cbind(arrowdf, catglab=mylabels)

# now create labels, arrows, etc.

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = catglab, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))

seaFiltPlot <- plot_ordination(sea, seaOrd, type = 'samples', color='site', title='sea')

seaFiltPlotArrow <- seaFiltPlot + 
  geom_segment(arrowmap, size = 1, data = arrowdf, color = "grey",  arrow = arrowhead, show_guide = FALSE) + 
  scale_color_manual(values=cols) +
  geom_text(labelmap, size = 2, data = arrowdf)
seaFiltPlotArrow


# see if SIMPER works

#metaFileSIMP <- subset(metaFile, site != 'NA')
seaDistSIMP <- merge(otu_table(sea), metaFileSIMP, by="row.names")

rownames(seaDistSIMP) <- seaDistSIMP[,1]
OTUsSIMP <- seaDistSIMP[,2:743]
factorsSIMP <- seaDistSIMP[,744:749]

spistSIMP <- simper(OTUsSIMP, factorsSIMP$site)

allTax["MED000008042", "Genus"]

# SIMPROF

seaSIMPROF <- simprof(OTUsSIMP, num.expected=10, num.simulated=9, method.cluster='average', method.distance='braycurtis', method.transform='squareroot', alpha=0.05, sample.orientation='row', silent=FALSE)

simprof.plot(seaSIMPROF, leafcolors=NA, plot=TRUE, fill=TRUE, leaflab="perpendicular", siglinetype=1)





