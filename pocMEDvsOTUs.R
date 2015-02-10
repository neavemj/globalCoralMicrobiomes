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

pocMED = read.table("poc.7801.matrixPercent", header=T)
rownames(pocMED) = pocMED[,1]
pocMED = pocMED[,2:length(pocMED)]

# Import normal taxonomy file from mothur

pocMEDTax = read.table('poc.7801.MED.nodeReps.nr_v119.knn.taxonomy', header=T, sep='\t')
rownames(pocMEDTax) = pocMEDTax[,1]
pocMEDTax = pocMEDTax[,3:9]
pocMEDTax = as.matrix(pocMEDTax)

# import meta data

metaFile = read.table('metaData2.MED', header=T, sep='\t')
rownames(metaFile) = metaFile[,1]
metaFile = metaFile[,2:7]

# 3% OTU and 1% OTU data for comparison

poc3OTUshared = read.table("poc.7801.0.03.0.03.pick.shared", header=T)
rownames(poc3OTUshared) = poc3OTUshared[,2]
poc3OTUshared = poc3OTUshared[,4:length(poc3OTUshared)]

poc1OTUshared = read.table("poc.7801.0.01.0.01.pick.shared", header=T)
rownames(poc1OTUshared) = poc1OTUshared[,2]
poc1OTUshared = poc1OTUshared[,4:length(poc1OTUshared)]

# Import taxonomy file from mothur

poc3OTUtax = read.table('poc.7801.unique.good.filter.precluster.an.0.03.cons.taxonomy', header=T, sep='\t')
rownames(poc3OTUtax) = poc3OTUtax[,1]
poc3OTUtax = poc3OTUtax[,3:9]
poc3OTUtax = as.matrix(poc3OTUtax)

poc1OTUtax = read.table('poc.7801.unique.good.filter.precluster.an.0.01.cons.taxonomy', header=T, sep='\t')
rownames(poc1OTUtax) = poc1OTUtax[,1]
poc1OTUtax = poc1OTUtax[,3:9]
poc1OTUtax = as.matrix(poc1OTUtax)

### Create phyloseq object

OTU = otu_table(pocMED, taxa_are_rows = FALSE)
TAX = tax_table(pocMEDTax) 
META = sample_data(metaFile)
pocPhylo = phyloseq(OTU, TAX, META)

OTUs3 = otu_table(poc3OTUshared, taxa_are_rows = FALSE)
OTUs1 = otu_table(poc1OTUshared, taxa_are_rows = FALSE)

TAX3 = tax_table(poc3OTUtax)
TAX1 = tax_table(poc1OTUtax)

poc3OTUphylo = phyloseq(OTUs3, TAX3, META)
poc1OTUphylo = phyloseq(OTUs1, TAX1, META)


# ordination comparing MED nodes and 3% OTUs

theme_set(theme_bw())
pocMEDord <- ordinate(pocPhylo, "NMDS", "bray")
plot_ordination(pocPhylo, pocMEDord, type = 'samples', color='site', title='MED nodes') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.23, 786 nodes

poc3OTUord <- ordinate(poc3OTUphylo, "NMDS", "bray")
plot_ordination(poc3OTUphylo, poc3OTUord, type = 'samples', color='site', title='3% OTUs') +
  geom_point(size=2) +
  scale_color_manual(values=cols)

# Best solution 0.23, 673 OTUs

poc1OTUord <- ordinate(poc1OTUphylo, "NMDS", "bray")
plot_ordination(poc1OTUphylo, poc1OTUord, type = 'samples', color='site', title='1% OTUs') +
  geom_point(size=2) +
  coord_flip() +
  scale_color_manual(values=cols)

# SAVE AS 700 x 532
# Best solution 0.23, 791 OTUs

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




