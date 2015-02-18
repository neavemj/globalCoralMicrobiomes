
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

metaFile = read.table('metaData2.MED', header=T, sep='\t')
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

# SAVE EPS 1500 x 800

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
plot_ordination(sea, seaOrd, type = 'samples', color='site', title='seawater')

# let's see if the chemistry is correlated with this
# first do PCA of just chemical data - can do this using R base packages
# import nutrient data

nutFile = read.table('Nutrients.txt', header=T, sep='\t')
rownames(nutFile) = nutFile[,1]
nutFile = nutFile[,2:(length(nutFile))]

PCA.res <- prcomp(nutFile, scale.=T)
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

site <- c("MaggieIs", "MaggieIs", "MaggieIs", "Ningaloo", "Ningaloo", "Ningaloo", "Ningaloo", "AmericanSamoa", "AmericanSamoa", "Maldives", "Maldives", "Maldives", "Maldives", "Maldives", "Maldives", "RedSea", "RedSea", "RedSea", "RedSea", "RedSea", "RedSea")
nutPCAggdata <- cbind(nutPCAggdata, site)

# create arrow info

arrowmap <- aes(xend = PC1, yend = PC2, x = 0, y = 0, alpha=0.5, shape = NULL, color = NULL, label = rownames(nutPCAarrows))
labelmap <- aes(x = PC1, y = PC2 + 0.04, shape = NULL, color = NULL, size=1.5, label = rownames(nutPCAarrows))
arrowhead = arrow(length = unit(0.02, "npc"))

nutPCA <- ggplot(nutPCAggdata, title='PCA of the nutrient data') +
  geom_point(aes(x=PC1, y=PC2, color=site), size=3) +
  scale_color_manual(values=c(cols)) 
nutPCA 

nutPCA + geom_segment(arrowmap, size = 0.5, data = nutPCAarrows, color = "black",  arrow = arrowhead, show_guide = FALSE) +
  geom_text(labelmap, size = 3, data = nutPCAarrows)



nutFile <- cbind(nutFile, site)

newMetaFile <- merge(metaFile, nutFile, by = "site", all.y = FALSE)



envfit(seaOrd, tmp$N.N)



# can also use Mothur to get Pearson correlations

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





