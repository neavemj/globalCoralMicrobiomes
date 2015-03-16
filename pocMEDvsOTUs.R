## Phyloseq and r analysis of the Minimum Entropy Decomposition results for poc ##
# 19.1.2015

# load required libraries

library("phyloseq")
library("ggplot2")
library("plyr")
library("vegan")
library("RColorBrewer")
setwd("./data")

# generate some colors to be consistent

#display.brewer.all()
#display.brewer.pal(n = 8, name = 'Dark2')
#brewer.pal(n = 8, name = "Dark2")

cols <- c("AmericanSamoa" = "#D95F02", "Indonesia" = "#A6761D", "MaggieIs" = "#666666", "Maldives" = "#E6AB02", "Micronesia" = "#66A61E", "Ningaloo" = "#7570B3", "RedSea" = "#E7298A")

# import normal percent matrix

allMED = read.table("all.7974.matrixPercent.txt", header=T)
rownames(allMED) = allMED[,1]
allMED = allMED[,2:length(allMED)]

# Import normal taxonomy file from mothur

allMEDTax = read.table('all.7974.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(allMEDTax) = allMEDTax[,1]
allMEDTax = allMEDTax[,3:9]
allMEDTax = as.matrix(allMEDTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

# 3% OTU and 1% OTU data for comparison

all3OTUshared = read.table("all.7974.0.03.pick.shared", header=T)
rownames(all3OTUshared) = all3OTUshared[,2]
all3OTUshared = all3OTUshared[,4:length(all3OTUshared)]

all1OTUshared = read.table("all.7974.0.01.pick.shared", header=T)
rownames(all1OTUshared) = all1OTUshared[,2]
all1OTUshared = all1OTUshared[,4:length(all1OTUshared)]

# Import taxonomy file from mothur

all3OTUtax = read.table('all.7974.0.03.taxonomy', header=T, sep='\t')
rownames(all3OTUtax) = all3OTUtax[,1]
all3OTUtax = all3OTUtax[,3:9]
all3OTUtax = as.matrix(all3OTUtax)

all1OTUtax = read.table('all.7974.0.01.taxonomy', header=T, sep='\t')
rownames(all1OTUtax) = all1OTUtax[,1]
all1OTUtax = all1OTUtax[,3:9]
all1OTUtax = as.matrix(all1OTUtax)

### Create phyloseq object

OTU = otu_table(allMED, taxa_are_rows = FALSE)
TAX = tax_table(allMEDTax) 
META = sample_data(metaFile)
allPhylo = phyloseq(OTU, TAX, META)

OTUs3 = otu_table(all3OTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(all1OTUshared, taxa_are_rows = FALSE)

TAX3 = tax_table(all3OTUtax)
TAX1 = tax_table(all1OTUtax)

all3OTUphylo = phyloseq(OTUs3, TAX3, META)
all1OTUphylo = phyloseq(OTUs1, TAX1, META)

pVerrPhylo <- subset_samples(allPhylo, species=='Pocillopora verrucosa')
pVerrPhylo = filter_taxa(pVerrPhylo, function(x) mean(x) > 0, TRUE)

pVerr3OTUphylo <- subset_samples(all3OTUphylo, species=='Pocillopora verrucosa')
pVerr3OTUphylo = filter_taxa(pVerr3OTUphylo, function(x) mean(x) > 0, TRUE)

pVerr1OTUphylo <- subset_samples(all1OTUphylo, species=='Pocillopora verrucosa')
pVerr1OTUphylo = filter_taxa(pVerr1OTUphylo, function(x) mean(x) > 0, TRUE)


# make mothur files into relative abundance

pVerr3OTUphyloRel = transform_sample_counts(pVerr3OTUphylo, function(x) x / sum(x) )
pVerr1OTUphyloRel = transform_sample_counts(pVerr1OTUphylo, function(x) x / sum(x) )

# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
pocMEDord <- ordinate(pVerrPhylo, "NMDS", "bray")
plot_ordination(pVerrPhylo, pocMEDord, type = 'samples', color='site', title='MED nodes') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.24, 762 nodes

poc3OTUord <- ordinate(pVerr3OTUphylo, "NMDS", "bray")
plot_ordination(pVerr3OTUphylo, poc3OTUord, type = 'samples', color='site', title='3% OTUs') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# Best solution 0.25, 680 OTUs

poc1OTUord <- ordinate(pVerr1OTUphylo, "NMDS", "bray")
plot_ordination(pVerr1OTUphylo, poc1OTUord, type = 'samples', color='site', title='1% OTUs') +
  geom_point(size=2) +
  coord_flip() +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.25, 833 OTUs

# do some bar plots to check if MED nodes resolve endozoic types better than OTUs

pocPhyloEndo = subset_taxa(pocPhylo, Genus=='Endozoicomonas')
poc3OTUphyloEndo = subset_taxa(poc3OTUphylo, Genus=='Endozoicomonas(100)')
poc1OTUphyloEndo = subset_taxa(poc1OTUphylo, Genus=='Endozoicomonas(100)')

pocPhyloEndoBar <- plot_bar(pocPhyloEndo, fill="Genus")
pocPhyloEndoBar + facet_wrap(~site, scales='free')

pocOTUphyloEndoBar <- plot_bar(poc3OTUphyloEndo, fill="Genus")
pocOTUphyloEndoBar + facet_wrap(~site, scales='free')

poc1OTUphyloEndoBar <- plot_bar(poc1OTUphyloEndo, fill="Genus")
poc1OTUphyloEndoBar + facet_wrap(~site, scales='free')

# add coloring for different Endozoicomonas OTUs

tax_table(pocPhyloEndo) <- cbind(tax_table(pocPhyloEndo), Strain=taxa_names(pocPhyloEndo))
tax_table(poc3OTUphyloEndo) <- cbind(tax_table(poc3OTUphyloEndo), Strain=taxa_names(poc3OTUphyloEndo))
tax_table(poc1OTUphyloEndo) <- cbind(tax_table(poc1OTUphyloEndo), Strain=taxa_names(poc1OTUphyloEndo))
myranks = c("Genus", "Strain")

mylabelsMED = apply(tax_table(pocPhyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels3OTU = apply(tax_table(poc3OTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")
mylabels1OTU = apply(tax_table(poc1OTUphyloEndo)[, myranks], 1, paste, sep="", collapse="_")

tax_table(pocPhyloEndo) <- cbind(tax_table(pocPhyloEndo), catglab=mylabelsMED)
tax_table(poc3OTUphyloEndo) <- cbind(tax_table(poc3OTUphyloEndo), catglab=mylabels3OTU)
tax_table(poc1OTUphyloEndo) <- cbind(tax_table(poc1OTUphyloEndo), catglab=mylabels1OTU)

# standardize the OTU shared files

poc3OTUphyloEndoFilt = filter_taxa(poc3OTUphyloEndo, function(x) mean(x) > 0.1, TRUE)
plot_bar(poc3OTUphyloEndoFilt, fill="catglab", title='3% OTUs') +
  facet_wrap(~site, scales='free')

poc1OTUphyloEndoFilt = filter_taxa(poc1OTUphyloEndo, function(x) mean(x) > 0.1, TRUE)
plot_bar(poc1OTUphyloEndoFilt, fill="catglab", title='1% OTUs') +
  facet_wrap(~site, scales='free')

pocPhyloEndoFilt = filter_taxa(pocPhyloEndo, function(x) mean(x) > 0.1, TRUE)
plot_bar(pocPhyloEndoFilt, fill="catglab", title='MED nodes') +
  facet_wrap(~site, scales='free')


# let's have a look at what happends to the most abundant Endo OTU

# 3% OTUs - OTU0001

poc3OTUphyloEndoRel = transform_sample_counts(poc3OTUphyloEndo, function(x) x / sum(x) )
poc3OTUphyloEndo1 = subset_taxa(poc3OTUphyloEndoRel, catglab=='Endozoicomonas(100)_Otu0001')
end1bar <- plot_bar(poc3OTUphyloEndo1, title='3% OTUs')

ggplot(end1bar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Maldives vs Red Sea

poc3OTUphyloEndo1Mal = subset_samples(poc3OTUphyloEndo1, site=='Maldives')
poc3OTUphyloEndo1RS = subset_samples(poc3OTUphyloEndo1, site=='RedSea')
poc3OTUphyloEndo1MalRS = merge_phyloseq(poc3OTUphyloEndo1Mal, poc3OTUphyloEndo1RS)

poc3OTUphyloEndo1MalRSBar <- plot_bar(poc3OTUphyloEndo1MalRS, title='3% OTUs')

ggplot(poc3OTUphyloEndo1MalRSBar$data, aes(x=site, y=Abundance, title='OTU1')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# 1% OTUs

poc1OTUphyloEndoRel = transform_sample_counts(poc1OTUphyloEndo, function(x) x / sum(x) )
pocOTUphyloEndo1_1 = subset_taxa(poc1OTUphyloEndoRel, catglab=="Endozoicomonas(100)_Otu00001")
pocOTUphyloEndo1_7 = subset_taxa(poc1OTUphyloEndoRel, catglab=="Endozoicomonas(100)_Otu00007")

pocOTUphyloEndo1_1 = merge_phyloseq(pocOTUphyloEndo1_1, pocOTUphyloEndo1_7)

end1bar_1 <- plot_bar(pocOTUphyloEndo1_1, title='1% OTUs', fill='catglab')

ggplot(end1bar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 


# just Maldives vs Red Sea

pocOTUphyloEndo1_1Mal = subset_samples(pocOTUphyloEndo1_1, site=='Maldives')
pocOTUphyloEndo1_1RS = subset_samples(pocOTUphyloEndo1_1, site=='RedSea')
pocOTUphyloEndo1_1MalRS = merge_phyloseq(pocOTUphyloEndo1_1Mal, pocOTUphyloEndo1_1RS)

pocOTUphyloEndo1_1MalRSBar <- plot_bar(pocOTUphyloEndo1_1MalRS, title='1% OTUs')

ggplot(pocOTUphyloEndo1_1MalRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab)) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

## SAVE AS 492 x 450

# OTU 1

pocOTUphyloEndo1_1Mal = subset_samples(pocOTUphyloEndo1_1, site=='Maldives')
pocOTUphyloEndo1_1RS = subset_samples(pocOTUphyloEndo1_1, site=='RedSea')
pocOTUphyloEndo1_1MalRS = merge_phyloseq(pocOTUphyloEndo1_1RS, pocOTUphyloEndo1_1Mal)

pocOTUphyloEndo1_1MalRSBar <- plot_bar(pocOTUphyloEndo1_1MalRS, title='OTU1')

ggplot(pocOTUphyloEndo1_1MalRSBar$data, aes(x=site, y=Abundance, title='OTU1')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450

# OTU 7

pocOTUphyloEndo1_7Mic = subset_samples(pocOTUphyloEndo1_7, site=='Maldives')
pocOTUphyloEndo1_7RS = subset_samples(pocOTUphyloEndo1_7, site=='RedSea')
pocOTUphyloEndo1_7MicRS = merge_phyloseq(pocOTUphyloEndo1_7RS, pocOTUphyloEndo1_7Mic)

pocOTUphyloEndo1_7MicRSBar <- plot_bar(pocOTUphyloEndo1_7MicRS, title='OTU7')

ggplot(pocOTUphyloEndo1_7MicRSBar$data, aes(x=site, y=Abundance, title='OTU7')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 


# MED nodes split from the original Otu0001 at 3%

pocMEDEndoRel = transform_sample_counts(pocPhyloEndo, function(x) x / sum(x) )

pocMEDEndoRel1_4919 = subset_taxa(pocMEDEndoRel, catglab=="Endozoicomonas_MED000004919")
pocMEDEndoRel1_6086 = subset_taxa(pocMEDEndoRel, catglab=="Endozoicomonas_MED000006086")
pocMEDEndoRel1_6080 = subset_taxa(pocMEDEndoRel, catglab=="Endozoicomonas_MED000006080")
pocMEDEndoRel1_5831 = subset_taxa(pocMEDEndoRel, catglab=="Endozoicomonas_MED000005831")
#pocMEDEndoRel1_4794 = subset_taxa(pocMEDEndoRel, catglab=="Endozoicomonas_MED000004794")

pocMEDEndoRel1 = merge_phyloseq(pocMEDEndoRel1_4919, pocMEDEndoRel1_6086, pocMEDEndoRel1_6080, pocMEDEndoRel1_5831)

endMEDbar_1 <- plot_bar(pocMEDEndoRel1, title='MED nodes', fill='catglab')
  
ggplot(endMEDbar_1$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

endMEDbar_6247 <- plot_bar(pocMEDEndoRel1_6247, title='MED nodes', fill='catglab')

ggplot(endMEDbar_6247$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# just Maldives vs Red Sea

pocMEDEndoRel1Mal = subset_samples(pocMEDEndoRel1, site=='Maldives')
pocMEDEndoRel1RS = subset_samples(pocMEDEndoRel1, site=='RedSea')
pocMEDEndoRel1MalRS = merge_phyloseq(pocMEDEndoRel1Mal, pocMEDEndoRel1RS)

pocMEDEndoRel1MalRSBar <- plot_bar(pocMEDEndoRel1MalRS, title='MED nodes')

ggplot(pocMEDEndoRel1MalRSBar$data, aes(x=site, y=Abundance)) +
  geom_boxplot(aes(fill=catglab), alpha=0.8) +
  geom_point(position=position_dodge(width=0.75), aes(group=catglab)) 

# MED6086

pocMEDEndoRel1_6086Mal = subset_samples(pocMEDEndoRel1_6086, site=='Maldives')
pocMEDEndoRel1_6086RS = subset_samples(pocMEDEndoRel1_6086, site=='RedSea')
pocMEDEndoRel1_6086MalRS = merge_phyloseq(pocMEDEndoRel1_6086RS, pocMEDEndoRel1_6086Mal)

pocMEDEndoRel1_6086MalRSBar <- plot_bar(pocMEDEndoRel1_6086MalRS, title='MED nodes')

ggplot(pocMEDEndoRel1_6086MalRSBar$data, aes(x=site, y=Abundance, title='MED6086')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols)

## SAVE AS 492 x 450

# MED6080

pocMEDEndoRel1_6080Mal = subset_samples(pocMEDEndoRel1_6080, site=='Maldives')
pocMEDEndoRel1_6080RS = subset_samples(pocMEDEndoRel1_6080, site=='RedSea')
pocMEDEndoRel1_6080MalRS = merge_phyloseq(pocMEDEndoRel1_6080RS, pocMEDEndoRel1_6080Mal)

pocMEDEndoRel1_6080MalRSBar <- plot_bar(pocMEDEndoRel1_6080MalRS, title='MED nodes')

ggplot(pocMEDEndoRel1_6080MalRSBar$data, aes(x=site, y=Abundance, title='MED6080')) +
  geom_boxplot(aes(fill = site)) +
  geom_point(position=position_jitter(width=0.2), aes(group=catglab, fill=site), size=5, pch=21, colour = 'black') +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) 

## SAVE AS 492 x 450




