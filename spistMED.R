
## Phyloseq and r analysis of the Minimum Entropy Decomposition results ##
# Stylophora pistillata samples extracted from all MED run
# 27.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library('ape')
library('labdsv')
library('grid')
library('RColorBrewer')
library("clustsig")
library("ellipse")
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

# import normal percent matrix

allShared = read.table("all.7974.matrixPercent.txt", header=T)
rownames(allShared) = allShared[,1]
allShared = allShared[,2:length(allShared)]

# Import normal taxonomy file from mothur

allTax = read.table('all.7974.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
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

# bar chart of s.pistillata phyla or class - the names parameter preserves order
# transforming to normal matix and as.factor keeps taxa stacked consistently

spist <- subset_samples(allPhylo, species=='Stylophora pistillata')

sample_data(spist)$names <- factor(sample_names(spist), levels=rownames(metaFile), ordered = TRUE)

spistFilt = filter_taxa(spist, function(x) mean(x) > 0.2, TRUE)

taxLevel <- "Class"

spistFiltGlom <- tax_glom(spistFilt, taxrank=taxLevel)
physeqdf <- psmelt(spistFiltGlom)

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


# SAVE PLOT: EPS 1500 x 600

#Top100OTUs = names(sort(taxa_sums(spist), TRUE)[1:100])
#spist100 = prune_taxa(Top50OTUs, spist)

# let's check what's happening with different Endozoicomonas OTUs

allPhyloEndo = subset_taxa(allPhylo, Genus=='Endozoicomonas')
spistEndo = subset_taxa(spist, Genus=='Endozoicomonas')

# add coloring for different Endozoicomonas OTUs

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), Strain=taxa_names(allPhyloEndo))
tax_table(spistEndo) <- cbind(tax_table(spistEndo), Strain=taxa_names(spistEndo))

myranks = c("Genus", "Strain")

mylabels = apply(tax_table(allPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabelsSpist = apply(tax_table(spistEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(allPhyloEndo) <- cbind(tax_table(allPhyloEndo), catglab=mylabels)
tax_table(spistEndo) <- cbind(tax_table(spistEndo), catglab=mylabelsSpist)

allPhyloEndoFilt = filter_taxa(allPhyloEndo, function(x) mean(x) > 0.1, TRUE)

sample_data(spistEndo)$names <- factor(sample_names(spistEndo), levels=rownames(metaFile), ordered = TRUE)

spistEndoFilt = filter_taxa(spistEndo, function(x) mean(x) > 0.2, TRUE)

theme_set(theme_bw())
plot_bar(spistEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_brewer(type='qual', palette = 'Dark2') +
  facet_grid(~site, scales='free', space='free_x')

# want to get different colors for endos in spist v poc
# first figure out how many diff endo types we have
coralEndos <- rbind(tax_table(spistEndoFilt)[,"catglab"], tax_table(pocEndoFilt)[,"catglab"])
unique(coralEndos)
length(unique(coralEndos))

# now get enough colors for each endo type
#display.brewer.all()
#display.brewer.pal(n = 12, name = 'Set3')
#brewer.pal(n = 12, name = 'Set3')
#display.brewer.pal(n = 8, name = 'Set1')
#brewer.pal(n = 8, name = 'Set1')

ggCols2 <- gg_color_hue(length(unique(coralEndos)))

cols2 <- c("Endozoicomonas_MED000010125" = "#F8766D", "Endozoicomonas_MED000009228"="#E38900", "Endozoicomonas_MED000009251"="#C49A00", "Endozoicomonas_MED000010206"="#A58AFF", "Endozoicomonas_MED000008196"="#53B400", "Endozoicomonas_MED000009620"="#00BC56", "Endozoicomonas_MED000007266"="#00C094", "Endozoicomonas_MED000005378"="#00BFC4", "Endozoicomonas_MED000010299"="#FF66A8", "Endozoicomonas_MED000009176"="#06A4FF", "Endozoicomonas_MED000009211"="#99A800", "Endozoicomonas_MED000009303"="#DF70F8", "Endozoicomonas_MED000001801"="#FB61D7", "Endozoicomonas_MED000009212"="#00B6EB")


theme_set(theme_bw())
plot_bar(spistEndoFilt, fill="catglab", x="names") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_manual(values=cols2) +
  facet_grid(~site, scales='free', space='free_x')

# SAVE PLOT: EPS 1500 x 600. Greater than 0.2%

# plot a heat map to show these differences a bit better

coralPhyloEndo = subset_samples(allPhyloEndo, species!='seawater')
spistPhyloEndo <- subset_samples(allPhyloEndo, species=='Stylophora pistillata')
pVerrPhyloEndo <- subset_samples(allPhyloEndo, species=='Pocillopora verrucosa')
spistPverrEndo <- merge_phyloseq(spistPhyloEndo, pVerrPhyloEndo)

spistPverrEndoFilt = filter_taxa(spistPverrEndo, function(x) mean(x) > 0.0, TRUE)

plot_heatmap(spistPverrEndoFilt, "RDA", "none", "species", "catglab", sample.order='species')
#plot_heatmap(spistPhyloEndo, "RDA", "none", "site", "catglab", sample.order='site')
#plot_heatmap(pVerrPhyloEndo, "RDA", "none", "site", "catglab", sample.order='site')

plot_heatmap(allPhylo, "RDA", "none", "species", "none")


# ordination for the spist samples

theme_set(theme_bw())
spistOrd <- ordinate(spist, "NMDS", "bray")
plot_ordination(spist, spistOrd, type = 'samples', color='site', title='spist')

# overlay the OTUs onto this

spistFilt = filter_taxa(spist, function(x) mean(x) > 0.1, TRUE)

#Top10OTUs = names(sort(taxa_sums(spistFilt), TRUE)[1:10])
#spistFilt10 = prune_taxa(Top10OTUs, spistFilt)

spistOrd <- ordinate(spist, "NMDS", "bray")
plot_ordination(spist, spistOrd, type = 'split', color='site', title='spist', label="Genus")

# calculate indicator species

spistDist <- vegdist(otu_table(spist), method="bray")
spistInd <- indspc(otu_table(spist), spistDist, numitr=1000)

scores.spist <- scores(spistSIMP, display='species')

spist.cor <- cor(otu_table(spist), method='pearson')

# sort by most important indicator species

spistIndVals <- spistInd$vals
spistIndValsSorted <- spistIndVals[order(-spistIndVals$indval),,drop=FALSE]

indval numocc  pval
MED000008130 0.5957302      2 0.003
MED000009492 0.5118680      4 0.001
MED000000375 0.5057752      2 0.008
MED000000581 0.4760191      2 0.011
MED000008713 0.4760191      2 0.011
MED000009640 0.4262102      5 0.001
MED000001212 0.4170521      2 0.021
MED000009470 0.4111436      4 0.001
MED000006552 0.4109851      3 0.001
MED000003841 0.3940892      9 0.001

# overlay these onto ordination

spistOrdTop10 <- spistOrd$species[c(rownames(spistIndValsSorted)[1:10]),]


spistMothurSpearman <- c("MED000006776",  "MED000009228",  "MED000005973",  "MED000009485",  "MED000005839",  "MED000007679",  "MED000009504",  "MED000009644", "MED000010328",  "MED000009696")


arrowmatrix <- data.frame(labels = spistOrdTop10)
#rownames(arrowmatrix) <- arrowmatrix[,1]
arrowdf <- data.frame(arrowmatrix)

# get taxonomic information from the original tax file

arrowdf <- data.frame(labels = allTax[rownames(arrowmatrix),"Genus"], arrowmatrix)
arrowdf <- cbind(arrowdf, spistOrd$species[spistMothurSpearman, ])

# also combine MED node with genus taxonomy for labelling

mylabels = apply(arrowdf[, c("labels", "labels.1")], 1, paste, sep="", collapse="_")
arrowdf <- cbind(arrowdf, catglab=mylabels)

# now create labels, arrows, etc.

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = arrowdf$labels)
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, label = arrowdf$labels, size=1.5)
arrowhead = arrow(length = unit(0.02, "npc"))

spistFiltPlot <- plot_ordination(spist, spistOrd, type = 'samples', color='site', title='spist') +
  scale_color_manual(values=cols)

spistFiltPlotArrow <- spistFiltPlot + 
  geom_segment(arrowmap, size = 1, data = arrowdf, color = "grey",  arrow = arrowhead, show_guide = FALSE) + 
  scale_color_manual(values=cols) +
  geom_text(labelmap, size = 2, data = arrowdf)

spistFiltPlotArrow


####################################################################
### fit chemistry to the coral microbiomes
####################################################################

waterQual <- c("temp", "salinity", "Domg", "pH")
nutrients <- c("PO4", "N.N", "silicate", "NO2", "NH4")
FCM <- c("prok", "syn", "peuk", "pe.peuk", "Hbact")

chemNoNA <- na.omit(metaFile[sample_names(spist),FCM])
spistNoNA <- prune_samples(rownames(chemNoNA), spist)

sample_names(spist)
sample_names(spistNoNA)

theme_set(theme_bw())
spistOrdNoNA <- ordinate(spistNoNA, "NMDS", "bray")
spistOrdNoNAPlot <- plot_ordination(spistNoNA, spistOrdNoNA, type = 'samples', color='site', title='spist') +
  geom_point(size=3) +
  scale_color_manual(values=c(cols)) 
spistOrdNoNAPlot

pointsNoNA <- spistOrdNoNA$points[rownames(chemNoNA),]

chemFit <- envfit(pointsNoNA, env = chemNoNA, na.rm=TRUE)

chemFit.scores <- as.data.frame(scores(chemFit, display= "vectors"))
chemFit.scores <- cbind(chemFit.scores, Species = rownames(chemFit.scores))

# create arrow info again

arrowmap <- aes(xend = MDS1, yend = MDS2, x = 0, y = 0, shape = NULL, color = NULL, label = rownames(chemFit.scores))
labelmap <- aes(x = MDS1, y = MDS2 + 0.04, shape = NULL, color = NULL, size=1.5, label = rownames(chemFit.scores))
arrowhead = arrow(length = unit(0.25, "cm"))

spistOrdNoNAPlot + 
  coord_fixed() +
  geom_segment(arrowmap, size = 0.5, data = chemFit.scores, color = "black",  arrow = arrowhead, show_guide = FALSE) +
  geom_text(labelmap, size = 3, data = chemFit.scores)

# SAVE AS 700 x 532


***VECTORS

MDS1     MDS2     r2 Pr(>r)    
temp      0.88707  0.46164 0.4174  0.001 ***
  salinity -0.15092 -0.98855 0.2369  0.005 ** 
  Domg     -0.71695 -0.69712 0.1044  0.057 .  
pH       -0.43162 -0.90206 0.3253  0.001 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.

***VECTORS

MDS1     MDS2     r2 Pr(>r)    
PO4      -0.40364  0.91492 0.1460  0.027 *  
  N.N       0.81490 -0.57960 0.0922  0.079 .  
silicate -0.75651  0.65398 0.5160  0.001 ***
  NO2      -0.55247  0.83353 0.4323  0.001 ***
  NH4       0.15821 -0.98740 0.1623  0.012 *  
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.

MDS1     MDS2     r2 Pr(>r)    
prok    -0.04407 -0.99903 0.2544  0.001 ***
  syn      0.51477  0.85733 0.1171  0.033 *  
  peuk     0.50215  0.86478 0.0906  0.065 .  
pe.peuk  0.95731  0.28907 0.0527  0.257    
Hbact   -0.73361  0.67957 0.0436  0.344    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
P values based on 999 permutations.


####################################################################
### SIMPROF analysis to check which samples fall into 'groups' without any *a priori* assumptions
####################################################################

Need to import the shared file containing just spist OTUs, then calcualte the simprof clusters based on the braycurtis metric. 

spistShared = otu_table(spist)
class(spistShared) <- "numeric"

spistSIMPROF <- simprof(spistShared, num.expected=10, num.simulated=9, method.cluster='average', method.distance='braycurtis', method.transform='identity', alpha=0.05, sample.orientation='row', silent=FALSE)

simprof.plot(spistSIMPROF, leafcolors=NA, plot=TRUE, fill=TRUE, leaflab="perpendicular", siglinetype=1)


# I'll try and overlay the significant clusters on top of the nMDS. 
# After calculating the clusters, make a data frame of the results and add to previous nMDS plot. Need to add these groups to the nMDS data.frame - I'll do a loop for this.

simprofCLUSTERS = data.frame()

for(j in 1:length(spistSIMPROF$significantclusters)){
  if(length(spistSIMPROF$significantclusters[[j]]) > 2){
    simprofCLUSTERS <- rbind(simprofCLUSTERS, cbind(j, spistSIMPROF$significantclusters[[j]]))
  }
}

rownames(simprofCLUSTERS) <- simprofCLUSTERS[,2]
colnames(simprofCLUSTERS) <- c("simprofCLUSTER", "group")
spistFiltDF <- as.data.frame(spistFiltPlotArrow$data)
metaSIMP <- merge(spistFiltDF, simprofCLUSTERS, by="row.names")

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

df_ell_saved <- (read.table("simpBrayIdent.txt"))
df_ell_saved$site <- factor(df_ell_saved$site)
nMDS.mean.saved <- read.table("simpBrayIdentMean.txt")

unique(df_ell$site)

wanted <- c("2", "3",  "5",  "8",  "10", "16", "21", "26", "27")

df_ell2 <- subset(df_ell, site==wanted)

spistFiltPlot +
  geom_path(data=df_ell, aes(x=x, y=y), color='black', size=0.5, linetype=2, show_guide=FALSE) +
  annotate("text", x=nMDS.mean$x, y=nMDS.mean$y, label=nMDS.mean$Group.1)

# SAVE AS 700 x 532

