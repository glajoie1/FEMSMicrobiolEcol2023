# Analyzing Field 2017 data

# Libraries
library(ggplot2)
library(ggvenn)
library(ape)
library(betapart)
library(devtools)
library(SpiecEasi)
library(biomformat)
library(tibble)
library(igraph)
library(bipartite)
library(car)
library(lme4)
library(phyloseq)
library(reshape2)
library(venn)
library(readr)

theme_set(theme_bw())

# Import objects
load("/home/lajoie/Garance/dada2_wk2.RData")

# Objects description
# ps.both.rar : epi and endo phyloseq object
# ps.epi.rar : epi only
# ps.endo.rar : endo only
# ps.pair.rar : epi and endo, but only including pairable samples (includes PRO and FRO)
# metadata_plotspp : metadata

##################
### Formatting ###
##################

# At least 3 samples/site/compartment dataset (n=69)

ps.acesac<-subset_samples(ps.both.rar, Species=='ACESAC')
# Drop JAC & PRO because less than 3 endophyte samples
ps.acesac<-subset_samples(ps.acesac, Site!='JAC')
ps.acesac<-subset_samples(ps.acesac, Site!='PRO')
ps.acesac<-prune_taxa(taxa_sums(ps.acesac)>0, ps.acesac)

# Paired dataset with 3 samples/site/compartment (n=42 samples)

ps.pair.acesac<-ps.pair.rar
# Remove PRO and FRO
ps.pair.acesac<-subset_samples(ps.pair.acesac, Site!='JAC')
ps.pair.acesac<-subset_samples(ps.pair.acesac, Site!='PRO')

# Generate 2 matching dataframes
ps.pair.join<-inner_join(subset_samples(ps.pair.acesac, Life_form=='Epiphyte')@sam_data, subset_samples(ps.pair.acesac, Life_form=='Endophyte')@sam_data, by=c('Site','Type','Plot'))
ps.pair.join<-do.call(cbind.data.frame, ps.pair.join)
ps.pair.sub<-ps.pair.join %>% group_by(Site) %>% arrange(Site,Type) %>% slice_head(n=3) #or replace arrange & slice_head by slice_sample(n = 3) # Random, so will not pick the same samples each time

# v3 as only epiphyte, v4 as only endophytes
v3<-subset_samples(ps.pair.acesac, Sample_ID%in%ps.pair.sub$Sample_ID.x)
v3<-prune_taxa(taxa_sums(v3)>0, v3)
v3@sam_data$pair<-paste(v3@sam_data$Site,v3@sam_data$Type,v3@sam_data$Plot, sep='-')
v4<-subset_samples(ps.pair.acesac, Sample_ID%in%ps.pair.sub$Sample_ID.y)
v4<-prune_taxa(taxa_sums(v4)>0, v4)
v4@sam_data$pair<-paste(v4@sam_data$Site,v4@sam_data$Type,v4@sam_data$Plot, sep='-')

#####################
## Alpha diversity ##
#####################

# All samples (n=83)

# Examine distribution of samples
table(ps.both.rar@sam_data$Life_form, ps.both.rar@sam_data$Site) # ACESAC: 32 samples of endophytes and 51 of epiphytes

# Calculate alpha diversity and add to metadata (richness + shannon)
sample_data(ps.both.rar)[,c('Richness','Shannon')] <- estimate_richness(ps.both.rar, measures=c('Observed','Shannon'))

# Differences in alpha-diversity across all samples (n=83)
adiv<-as.data.frame(ps.both.rar@sam_data)
leveneTest(adiv$Richness~adiv$Life_form) # Variance not homogeneous, take non-parametric test
wilcox.test(adiv$Richness) # Sig p<0.001
wilcox.test(adiv$Shannon) # Sig p<0.001

# Plot
adiv.plot<-melt(adiv, measure.vars=c('Richness','Shannon'))
ggplot(adiv.plot, aes(x=Life_form, y=value))+
  geom_boxplot() +
  facet_wrap(~variable, scale='free')+
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=12), strip.text = element_text(size=12))

# On paired samples (n=56)

# Effect of life form and site on richness

# Calculate alpha diversity and add to metadata (richness + shannon)
sample_data(ps.pair.acesac)[,c('Richness','Shannon')] <- estimate_richness(ps.pair.acesac, measures=c('Observed','Shannon'))

adiv.pair<-as.data.frame(ps.pair.acesac@sam_data)
leveneTest(adiv.pair$Richness~adiv.pair$Life_form*adiv.pair$Site) # variances are homogeneous
leveneTest(adiv.pair$Shannon~adiv.pair$Life_form*adiv.pair$Site) # variances are homogeneous

# Effect of Life_form & Site

# Richness
rich.t<-aov(adiv.pair$Richness~adiv.pair$Life_form*adiv.pair$Site)
#rich.t<-aov(meta.new$Richness~meta.new$Life_form*meta.new$MAT)
summary(rich.t)
plot(rich.t)

# Shannon
shan.t<-aov(adiv.pair$Shannon~adiv.pair$Life_form*adiv.pair$Site)
# shan.t<-aov(meta.new$Shannon~meta.new$Life_form*meta.new$MAT)
summary(shan.t)
plot(shan.t)

############
# Taxonomy #
############

#### Taxonomic composition
# All 83 samples
p1<-subset_samples(ps.both.rar, Life_form=='Epiphyte')
p2<-subset_samples(ps.both.rar, Life_form=='Endophyte')

## Epiphytes

# Append OTU table with taxo table
otu.tab<-as.data.frame(t(p1@otu_table@.Data))
# Sum of each ASV freq across samples
otu.sum<-as.data.frame(rowSums(otu.tab))
# Append (average across samples)
taxo16<-as.data.frame(p1@tax_table@.Data)
taxo16[is.na(taxo16)]<-'Z-Unidentified'
#taxo16$phylum[which(is.na(taxo16$phylum)==T)]<-'Unidentified phylum'
otu.taxo16<-merge(otu.sum,taxo16, by='row.names')
rownames(otu.taxo16)<-otu.taxo16$Row.names
otu.taxo16<-otu.taxo16[,-1]
colnames(otu.taxo16)[1]<-'Freq'
# Aggregate by phylum
phylo_sum<-aggregate(Freq~class, otu.taxo16, sum)
# Add unidentified bacteria
phylo_sum$Freq<-as.numeric(phylo_sum$Freq)
# Relative abundance
phylo_sum$Freq<-phylo_sum$Freq/sum(phylo_sum$Freq)
phylo_sum$Life_form<-'Epiphyte'

## Endophytes

# Append OTU table with taxo table
otu.tab<-as.data.frame(t(p2@otu_table@.Data))
# Sum of each ASV freq across samples
otu.sum<-as.data.frame(rowSums(otu.tab))
# Append (average across samples)
taxo16<-as.data.frame(p2@tax_table@.Data)
taxo16[is.na(taxo16)]<-'Z-Unidentified'
#taxo16$phylum[which(is.na(taxo16$phylum)==T)]<-'Unidentified phylum'
otu.taxo16<-merge(otu.sum,taxo16, by='row.names')
rownames(otu.taxo16)<-otu.taxo16$Row.names
otu.taxo16<-otu.taxo16[,-1]
colnames(otu.taxo16)[1]<-'Freq'
# Aggregate by phylum
phylo_sum2<-aggregate(Freq~class, otu.taxo16, sum)
# Add unidentified bacteria
phylo_sum2$Freq<-as.numeric(phylo_sum2$Freq)
# Relative abundance
phylo_sum2$Freq<-phylo_sum2$Freq/sum(phylo_sum2$Freq)
phylo_sum2$Life_form<-'Endophyte'

# Merge
phylo_sum_both<-rbind(phylo_sum, phylo_sum2)

# Plot -> starts to vary more at the level of order
# Need to make sure that 'non-assigned' show up
ggplot(phylo_sum_both, aes(x=Life_form, y=Freq, fill=class))+
  geom_bar(color='black', stat='identity',position='stack')+
  theme_bw()+
  guides(fill = guide_legend(ncol = 1))+
  labs(x='Life form', y='Relative abundance')


# Per site, per compartment taxo

# Paired only without JAC, PRO ; n=56 samples
p1<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte')
p2<-subset_samples(ps.pair.acesac, Life_form=='Endophyte')

# 3 paired per site (n=42 samples)
# p1=v3
# p2=v4

# Append OTU table with taxo table
otu.tab<-melt(p1@otu_table@.Data)
otu.tab<-merge(otu.tab, metadata_plotspp[,c('Sample_ID','Site')], by.x='Var1', by.y='Sample_ID', all.x=T)
# Sum per site
otu.site<-aggregate(value~Var2+Site,otu.tab, sum)

# Append (average across samples)
taxo16<-as.data.frame(p1@tax_table@.Data)
taxo16[is.na(taxo16)]<-'Z-Unidentified'
otu.taxo16<-merge(otu.site,taxo16, by.x='Var2', by.y='row.names', all.x=T)
colnames(otu.taxo16)[3]<-'Freq'

# Aggregate by phylum
phylo_sum<-aggregate(Freq~class+Site, otu.taxo16, sum)
# Add unidentified bacteria
phylo_sum$Freq<-as.numeric(phylo_sum$Freq)
site.sum<-aggregate(Freq~Site, phylo_sum, sum)
colnames(site.sum)[2]<-'Total_freq'
phylo_sum<-merge(phylo_sum, site.sum, by='Site', all.x=T)
# Relative abundance
phylo_sum$Freq<-phylo_sum$Freq/phylo_sum$Total_freq #(3 samples at 4000 seq per site)
phylo_sum$Life_form<-'Epiphyte'

## Endophytes

# Append OTU table with taxo table
otu.tab<-melt(p2@otu_table@.Data)
otu.tab<-merge(otu.tab, metadata_plotspp[,c('Sample_ID','Site')], by.x='Var1', by.y='Sample_ID', all.x=T)
# Sum per site
otu.site<-aggregate(value~Var2+Site,otu.tab, sum)

# Append (average across samples)
taxo16<-as.data.frame(p2@tax_table@.Data)
taxo16[is.na(taxo16)]<-'Z-Unidentified'
otu.taxo16<-merge(otu.site,taxo16, by.x='Var2', by.y='row.names', all.x=T)
colnames(otu.taxo16)[3]<-'Freq'

# Aggregate by phylum
phylo_sum2<-aggregate(Freq~class+Site, otu.taxo16, sum)
# Add unidentified bacteria
phylo_sum2$Freq<-as.numeric(phylo_sum2$Freq)
site.sum<-aggregate(Freq~Site, phylo_sum2, sum)
colnames(site.sum)[2]<-'Total_freq'
phylo_sum2<-merge(phylo_sum2, site.sum, by='Site', all.x=T)
# Relative abundance
phylo_sum2$Freq<-phylo_sum2$Freq/phylo_sum2$Total_freq #(3 samples at 4000 seq per site)
phylo_sum2$Life_form<-'Endophyte'

# Merge
phylo_sum_both<-rbind(phylo_sum, phylo_sum2)

phylo_sum_both$Site<- factor(phylo_sum_both$Site, levels=c('VAL','SMS','VER','MEG','GAT','MSH','FRO'))

# Plot -> starts to vary more at the level of order
# Need to make sure that 'non-assigned' show up
ggplot(phylo_sum_both, aes(x=Life_form, y=Freq, fill=class))+
  geom_bar(color='black', stat='identity',position='stack')+
  theme_bw()+
  facet_wrap(~Site, nrow=2)+
  guides(fill = guide_legend(ncol = 1))+
  labs(x='Life form', y='Relative abundance')


###################
## Venn diagrams ##
###################

# How many ASVs are shared among endophytes and epiphytes
# Do this using only paired samples
# Venn diagram
v1<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte') # & Site=='VAL'
v2<-subset_samples(ps.pair.acesac, Life_form=='Endophyte')

v12<-list(Endophytes=names(which(colSums(v2@otu_table)>0)),Epiphytes=names(which(colSums(v1@otu_table)>0)))
ggvenn(v12,stroke_size = 0.3, set_name_size = 5)
#venn(v12)

table(ps.pair.acesac@sam_data$Life_form,ps.pair.acesac@sam_data$Site)

####
# Determine the number of OTUs that are part of epi, shared, or total for each pair, and generate statistis

ps.pair.acesac@sam_data$pair<-paste(ps.pair.acesac@sam_data$Site,ps.pair.acesac@sam_data$Type,ps.pair.acesac@sam_data$Plot, sep='-')

pair.list<-NULL
for (i in unique(ps.pair.acesac@sam_data$pair)){
  vep<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte' & ps.pair.acesac@sam_data$pair==i) 
  ven<-subset_samples(ps.pair.acesac, Life_form=='Endophyte' & ps.pair.acesac@sam_data$pair==i) 
  
  v12<-list(Endophytes=names(which(colSums(ven@otu_table)>0)),Epiphytes=names(which(colSums(vep@otu_table)>0)))
  
  # print(ggvenn(v12,stroke_size = 0.3, set_name_size = 5)+
  #        ggtitle(i))
  pair.list<-rbind(pair.list,c(length(setdiff(v12$Endophytes,v12$Epiphytes)), length(intersect(v12$Epiphytes, v12$Endophytes)), length(setdiff(v12$Epiphytes,v12$Endophytes)), length(union(v12$Epiphytes,v12$Endophytes))))
  
  #Sys.sleep(3)
}
rownames(pair.list)<-unique(ps.pair.acesac@sam_data$pair)
colnames(pair.list)<-c('endo','shared','epi','total')
pair.list<-cbind(pair.list, colsplit(rownames(pair.list),'-', names=c('Site','Type','Plot')))
pair.list$pair<-rownames(pair.list)
pair.list.m<-melt(pair.list, measure.vars=c('endo','shared','epi','total'))

ggplot(pair.list.m, aes(y=value, x=variable))+ # , color=Site
  geom_boxplot()+
 # geom_point(position=position_dodge(width=0.75),aes(group=variable))+
  scale_y_continuous(breaks = round(seq(min(pair.list.m$value), max(pair.list.m$value), by = 50)))+
  theme_bw()

# Some stats
mean(pair.list.m[which(pair.list.m$variable=='shared'),'value'])
mean(pair.list.m[which(pair.list.m$variable=='endo'),'value'])
pair.list$prop.shared<-pair.list$shared/pair.list$total
mean(pair.list$prop.shared)

# Differences among mean values - repeated measures
pair.list.m$pair<-as.factor(pair.list.m$pair)
leveneTest(pair.list.m$value~pair.list.m$variable) # Variance not homogeneous, take non-parametric test
anova_test(pair.list.m, value~variable) # Variance not homogeneous, take non-parametric test
friedman.test(pair.list.m$value, pair.list.m$variable, pair.list.m$pair) # Sig p<0.001
# Order variables for pairing in wilcox test
#https://www.datanovia.com/en/lessons/friedman-test-in-r/
pair.list.m.ord<-pair.list.m[order(pair.list.m$variable,pair.list.m$pair),c(4,5,6)]
pair.list.m.ord$value<-as.numeric(pair.list.m.ord$value)
pair.list.m.ord<-as_tibble(pair.list.m.ord)
#pair.list.m.ord<-pair.list.m.ord[-which(pair.list.m.ord$variable%in%c('total')),]

pair.list.m.ord %>% friedman_test(value~variable|pair)
pairwise.wilcox.test(x=pair.list.m.ord$value,g=pair.list.m.ord$variable, paired=TRUE, exact=F)



# Based on abundance
pair.abund<-NULL
for (i in unique(ps.pair.acesac@sam_data$pair)){
  vep<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte' & ps.pair.acesac@sam_data$pair==i) 
  ven<-subset_samples(ps.pair.acesac, Life_form=='Endophyte' & ps.pair.acesac@sam_data$pair==i) 
  
  v12<-list(Endophytes=names(which(colSums(ven@otu_table)>0)),Epiphytes=names(which(colSums(vep@otu_table)>0)))
  
  gg<-intersect(v12$Endophytes,v12$Epiphytes)
  if(length(gg)>0){
    # OPTION A
  gepi<-sum(vep@otu_table[,which(colnames(vep@otu_table)%in%gg)])/4000
  gedo<-sum(ven@otu_table[,which(colnames(ven@otu_table)%in%gg)])/4000

  pair.abund<-rbind(pair.abund,c(gepi,gedo))


  }

}

rownames(pair.abund)<-unique(ps.pair.acesac@sam_data$pair)
colnames(pair.abund)<-c('epi','endo')


### Venn diagram among sites

# Epiphytes
par(mfrow = c(2, 4))

# Generate list per site
epilist<-list()
for (i in 1:7){
  sub<-subset_samples(v3, Site%in%unique(Site)[i])
  epilist[[i]]<-names(which(colSums(sub@otu_table)>0))
}

names(epilist)<-unique(v3@sam_data$Site)
venn(epilist, ilcs=0.8)


uni.plot<-NULL
# Those that are shared among plots within a site are not unique to that site
for (i in 1:7){
  
sub<-subset_samples(v3, Site%in%unique(Site)[i])

i1<-names(which(colSums(sub@otu_table[1,])>0))
i2<-names(which(colSums(sub@otu_table[2,])>0))
i3<-names(which(colSums(sub@otu_table[3,])>0))

ffr<-setdiff(unlist(epilist[[i]]),unlist(epilist[-i])) 
list.end.s<-list(i1,i2,i3) 

# Prop unique to plots
uni.plot<-c(uni.plot, length(c(Reduce(setdiff, list(i1,i2,i3)),Reduce(setdiff, list(i2,i3,i1)),Reduce(setdiff, list(i3,i1,i2))))/length(Reduce(union, list(i1,i2,i3))))


}

mean(uni.plot) # 65%



list.end.s<-list(i1,i2,i3,ffr)
venn(list.end.s)
int.site<-Reduce(intersect, list(i1,i2,i3))


# Check which epiphytes are common to all sites
Reduce(intersect, epilist) # 115 ASVs
# Who are they
taxcore<-data.frame(v3@tax_table[Reduce(intersect, epilist)])
core.epi<-table(taxcore$family) # , useNA='ifany'
# Are they overrepresented? Meh, all around the same proportion...
gen.epi<-table(v3@tax_table[,6])
gen.epi<-gen.epi[which(names(gen.epi)%in%names(table(taxcore$genus, useNA='ifany')))]
core.epi/gen.epi
sum(core.epi)/sum(gen.epi)

# Total number of ASVs
table(taxa_sums(v3)>0) # 1736 ASVs
# Prop of shared ASV
length(Reduce(intersect, epilist))/table(taxa_sums(v3)>0)[[1]] # ~ 6.6%
# How much rel. abundance
sum(v3@otu_table[,which(colnames(v3@otu_table)%in%Reduce(intersect, epilist))]) /sum(v3@otu_table) # ~65% of total number of seq across sites


# # of unique ASVs for each site
uni.asv<-NULL
for (i in 1:7){
  uni.asv<-c(uni.asv, length(setdiff(unlist(epilist[i]),unlist(epilist[-i]))))
}
names(uni.asv)<-names(epilist)
uni.asv
# Sum & Average
sum(uni.asv)
mean(uni.asv) # ~135 ASVs in average are unique to sites

# Endophytes

# Generate list per site
par(mfrow = c(2, 4))

endlist<-list()
for (i in 1:7){
  sub<-subset_samples(v4, Site%in%unique(Site)[i])
  endlist[[i]]<-names(which(colSums(sub@otu_table)>0))
}

names(endlist)<-unique(v4@sam_data$Site)
venn(endlist, ilcs=0.8)

lapply(endlist, FUN=length)

# Those that are shared among plots within a site are not unique to that site
uni.plot<-NULL

for (i in 1:7){
sub<-subset_samples(v4, Site%in%unique(Site)[i])

i1<-names(which(colSums(sub@otu_table[1,])>0))
i2<-names(which(colSums(sub@otu_table[2,])>0))
i3<-names(which(colSums(sub@otu_table[3,])>0))


ffr<-setdiff(unlist(endlist[[i]]),unlist(endlist[-i])) # Unique to MEG
list.end.s<-list(i1,i2,i3) # All endophytes at MEG
#venn(list.end.s, ilcs=1.5)

# Prop unique to plots
uni.plot<-c(uni.plot, length(c(Reduce(setdiff, list(i1,i2,i3)),Reduce(setdiff, list(i2,i3,i1)),Reduce(setdiff, list(i3,i1,i2))))/length(Reduce(union, list(i1,i2,i3))))

}

names(uni.plot)<-unique(v4@sam_data$Site)

mean(uni.plot) # 90%

# Total number of ASVs
table(taxa_sums(v4)>0) # 515 ASVs
# Prop of shared ASV
length(Reduce(intersect, endlist))/table(taxa_sums(v4)>0)[[1]] # ~ 0.3%
# How much rel. abundance
sum(v4@otu_table[,which(colnames(v4@otu_table)%in%Reduce(intersect, endlist))]) /sum(v4@otu_table) # ~8.6% of total number of seq across sites


list.end.s<-list(i1,i2,i3,ffr)
venn(list.end.s)
int.site<-Reduce(intersect, list(i1,i2,i3))


# Unique to plots / total
length(c(Reduce(setdiff, list(i1,i2,i3)),Reduce(setdiff, list(i2,i3,i1)),Reduce(setdiff, list(i3,i1,i2))))/length(Reduce(union, list(i1,i2,i3)))


unlist(lapply(endlist, FUN=length))

# Check which are common to all sites
Reduce(intersect, endlist) # 2 ASVs
# Who are they
v4@tax_table[Reduce(intersect, endlist)] # 2 beijerinckiaceae

# Total number of ASVs
table(taxa_sums(v4)>0) # 526 ASVs
# Prop
length(Reduce(intersect, endlist))/table(taxa_sums(v4)>0)[[1]] # ~ 0.4%

# # of unique ASVs for each site
uni.asv2<-NULL
for (i in 1:7){
  uni.asv2<-c(uni.asv2, length(setdiff(unlist(endlist[i]),unlist(endlist[-i]))))
}
names(uni.asv2)<-names(endlist)
# Sum & Average
sum(uni.asv2)
mean(uni.asv2) # 56 ASVs

####################
## Beta-diversity ##
####################

#### Ordination

# Across all samples
# Bray-Curtis distance on OTU matrix
mcc<-vegdist(otu_table(ps.acesac), method='bray')
# PCoA
mcc.pcoa <-cmdscale(mcc, k=nrow(otu_table(ps.acesac))-1, eig=T)
mcc.pcoa
# Eigenvalues
mcc.pcoa$eig/sum(mcc.pcoa$eig) # 16.0% first axis, 7.3% second axis, 4.8% third axis

groups=ps.acesac@sam_data$Life_form
color.points<-c('black','white')[as.numeric(as.factor(groups))]

# MAT
efit.mat<-envfit(mcc.pcoa, ps.acesac@sam_data[,c('MAT','MAP'), drop=F], choices=c(1:2)) # Only MAT sig

# Plot
ord<-ordiplot(mcc.pcoa, choices=c(1,2),type='none') # , type='text', display='site'
points(ord, 'sites', cex=1, pch=21, bg=color.points)
ordiellipse(mcc.pcoa, choices=c(1,2), groups, display = "sites", kind = "se", lty=c(3,1), conf = 0.95, alpha = 0.05, lwd = 1.5, label=TRUE)
plot(efit.mat, lty=3, col='black', cex=0.8)


#Plot By Site

groups2=ps.acesac@sam_data$Site
color.points<-c('blue','red','yellow','green','orange','grey','pink')[as.numeric(as.factor(groups2))]

ord<-ordiplot(mcc.pcoa, choices=c(2,3),type='none') # , type='text', display='site'
points(ord, 'sites', cex=0.9, pch=19, col=color.points)
ordispider(mcc.pcoa, choices=c(2,3), groups2, display = "sites", kind = "se", lty=c(1), conf = 0.95, alpha = 0.05, lwd = 1.5, label=TRUE, col=c('blue','red','yellow','green','orange','grey','pink'))
plot(efit.mat, lty=3, col='black', cex=0.8)

# PERMANOVA for effect of life form and site
di<-betadisper(mcc, ps.acesac@sam_data$Life_form)
anova(di) # Variance not homogeneous among life forms
permutest(di, perm=999)
perm.lifesite<-adonis2(ps.pair.acesac@otu_table~Life_form*Site, data.frame(ps.pair.acesac@sam_data))
perm.lifesite

# Is effect of site greater for epiphyte or endophyte
fgh<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte')
di<-betadisper(vegdist(fgh@otu_table, 'bray'), fgh@sam_data$Site)
anova(di)
adonis2(fgh@otu_table~Site, data.frame(fgh@sam_data))

fgh<-subset_samples(ps.pair.acesac, Life_form=='Endophyte')
di<-betadisper(vegdist(fgh@otu_table, 'bray'), fgh@sam_data$Site)
anova(di)
adonis2(fgh@otu_table~Site, data.frame(fgh@sam_data))


######################
## Deseq analysis ##
####################

# https://joey711.github.io/phyloseq-extensions/DESeq2.html

library("DESeq2")
enrich.life<-phyloseq_to_deseq2(ps.acesac, ~Life_form)
enrich.life<-DESeq(enrich.life, test='Wald', fitType = 'parametric')

# Problems with 0s
# https://github.com/joey711/phyloseq/issues/1556
#https://github.com/joey711/phyloseq/issues/642
#https://github.com/joey711/phyloseq/issues/445
# So first apply filtering to keep only taxa that have >20 counts
ps.deseq<-ps.acesac
ps.deseq<-filter_taxa(ps.deseq, function(x) sum(x)>=50, TRUE)
colSums(ps.deseq@otu_table)

enrich.life<-phyloseq_to_deseq2(ps.deseq, ~Life_form)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans = apply(counts(enrich.life), 1, gm_mean)

dds1 = estimateSizeFactors(enrich.life, geoMeans=geoMeans)
dds1<-estimateDispersions(dds1)
dds1 <- nbinomWaldTest(dds1)

res = results(dds1, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
# Log fold change positive (enriched in epiphytes)
# Log fold change negative (enriched in endophytes)

# Plotting them

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Class order
x = tapply(sigtab$log2FoldChange, sigtab$class, function(x) max(x))
x = sort(x, TRUE)
sigtab$class = factor(as.character(sigtab$class), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=class)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Very few are enriched in the endophytes
# No clear separation in the genera where there are enriched ASVs among compartments + no specific phylum that are overrepresented (a bit of everything)
table(sigtab$class)



#### Nestedness vs turnover : Endo vs epi

# Performed on paired samples

# EPI table
epi.tab<-ps.pair.acesac@otu_table@.Data
epi.tab<-epi.tab[match(ps.pair.sub$Sample_ID.x,rownames(epi.tab)),]
epi.tab.pa<-decostand(epi.tab[match(ps.pair.sub$Sample_ID.x,rownames(epi.tab)),], 'pa')
rownames(epi.tab.pa)<-paste0('EE',seq(1:nrow(epi.tab.pa)))

# Endo table
endo.tab<-ps.pair.acesac@otu_table@.Data
endo.tab<-endo.tab[match(ps.pair.sub$Sample_ID.y,rownames(endo.tab)),]
endo.tab.pa<-decostand(endo.tab[match(ps.pair.sub$Sample_ID.y,rownames(endo.tab)),], 'pa')
rownames(endo.tab.pa)<-paste0('EE',seq(1:nrow(epi.tab.pa)))

# Beta partitioning analysis
###### Paired: presence absence only
#beta.sim -> turnover component
#beta.sne -> nestedness component
#beta.sor -> overall 
betapart.both<-beta.temp(epi.tab.pa, endo.tab.pa, 'sorensen')
betapart.both$epi<-rownames(epi.tab)
betapart.both$endo<-rownames(endo.tab)
# Join metadata
betapart.both<-left_join(betapart.both, metadata_plotspp, by=c('epi'='Sample_ID'))
betapart.both$Site<-as.factor(betapart.both$Site)
betapart.both$beta.prop.sim<-betapart.both$beta.sim/betapart.both$beta.sor

hist(betapart.both$beta.sim/betapart.both$beta.sor)
summary(betapart.both$beta.sim/betapart.both$beta.sor)
sd(betapart.both$beta.sim/betapart.both$beta.sor)
hist(betapart.both$beta.sor)


###############################
## Ordination by compartment ##
###############################

# Just epiphytes

# Community matrix
otu.epi<-subset_samples(ps.acesac, Life_form=='Epiphyte')
otu.epi.otu<-otu.epi@otu_table@.Data
otu.epi.otu<-decostand(otu.epi.otu, 'hellinger')

# Geographic positions
library(geodist)
geo<-metadata_plotspp[which(metadata_plotspp$Sample_ID%in%rownames(otu.epi.otu)),]
geode<-data.frame(lon=geo$Coord.x, lat=geo$Coord.y)
geode.di<-geodist(geode, measure='geodesic')
rownames(geode.di)<-colnames(geode.di)<-geo$Sample_ID
geode.di<-geode.di[match(rownames(otu.epi.otu),rownames(geode.di)),match(rownames(otu.epi.otu),colnames(geode.di))]
geode.di<-as.dist(geode.di, diag=T)

otu.epi.dist<-vegdist(otu.epi@otu_table@.Data, 'bray', diag=T)

mantel.test(as.matrix(otu.epi.dist),as.matrix(geode.di))


# Plant traits
otu.epi.trt<-otu.epi@sam_data
otu.epi.trt<-otu.epi.trt[,which(colnames(otu.epi.trt)%in%c('SLA','P','Ca','Wood.dens','MAT','MAP'))]
otu.epi.trt<-as.data.frame(as.matrix(otu.epi.trt))

# Standardize traits
# Fill NA values with the mean of the site-level measurements
otu.epi.trt$Wood.dens[is.na(otu.epi.trt$Wood.dens)==T]<-c(0.625,0.619, 0.606)
otu.epi.trt.std<-decostand(as.data.frame(otu.epi.trt), 'standardize')
plot(match(rownames(otu.epi.trt.std), rownames(otu.epi.otu)))
# RDA with traits as constraining variables
epi.rda.trt<-rda(otu.epi.otu, otu.epi.trt.std[,1:4], otu.epi.trt.std[,5:6])
epi.rda.trt<-rda(otu.epi.otu~., otu.epi.trt.std[,1:6])
# Site colors
col.site<-as.numeric(as.factor(metadata_plotspp$Site[match(labels(epi.rda.trt$CCA$u)[[1]],metadata_plotspp$Sample_ID)]))
shp.site<-c(3,8,21,22,23,24,25)[col.site]
# Plot with site as color
plot(epi.rda.trt, type='n',)
points(epi.rda.trt, 'sites', pch=shp.site, cex=1.4, bg=col.site)
text(epi.rda.trt, 'bp', col='red', cex=1.0, arrow.mul=1.0)


summary(epi.rda.trt)

ordiR2step(rda(otu.epi.otu~1, data=otu.epi.trt.std[,1:6],na.action=na.exclude), 
           scope= formula(epi.rda.trt), direction= "forward", R2scope=TRUE, pstep=1000)

adonis2(otu.epi.otu~ MAT + SLA + Ca, data=otu.epi.trt.std, method='bray', permutations=999) # Changing the order changes sig.
adonis2(otu.epi.otu~., otu.epi.trt.std[,1:6], method='bray', permutations=999)

# Variation partitioning
var.epi<-varpart(otu.epi.otu, otu.epi.trt.std[,1:4], otu.epi.trt.std[,5:6]) #X1 is plant traits, X2 is environmental variables
plot(var.epi)

# Just endophytes -> Use rda with plant traits as constraining variables
# Community matrix
otu.endo<-subset_samples(ps.acesac, Life_form=='Endophyte')
otu.endo.otu<-otu.endo@otu_table@.Data
otu.endo.otu<-decostand(otu.endo.otu, 'hellinger')

library(geodist)
geo<-metadata_plotspp[which(metadata_plotspp$Sample_ID%in%rownames(otu.endo.otu)),]
geode<-data.frame(lon=geo$Coord.x, lat=geo$Coord.y)
geode.di<-geodist(geode, measure='geodesic')
rownames(geode.di)<-colnames(geode.di)<-geo$Sample_ID
geode.di<-geode.di[match(rownames(otu.endo.otu),rownames(geode.di)),match(rownames(otu.endo.otu),colnames(geode.di))]
geode.di<-as.dist(geode.di, diag=T)

otu.endo.dist<-vegdist(otu.endo@otu_table@.Data, 'bray', diag=T)

mantel.test(as.matrix(otu.endo.dist),as.matrix(geode.di))



# Plant traits
otu.endo.trt<-otu.endo@sam_data
otu.endo.trt<-otu.endo.trt[,which(colnames(otu.endo.trt)%in%c('SLA','P','Ca','Wood.dens','MAT','MAP'))]
otu.endo.trt<-as.data.frame(as.matrix(otu.endo.trt))
# Standardize traits
# Fill NA values with the mean of the site-level measurements
otu.endo.trt$Wood.dens[is.na(otu.endo.trt$Wood.dens)==T]<-c(0.625,0.606)
otu.endo.trt.std<-decostand(as.data.frame(otu.endo.trt), 'standardize')
plot(match(rownames(otu.endo.trt.std), rownames(otu.endo.otu)))
# RDA with traits as constraining variables
endo.rda.trt<-rda(otu.endo.otu, otu.endo.trt.std[,1:4], otu.endo.trt.std[,5:6])
endo.rda.trt<-rda(otu.endo.otu~., otu.endo.trt.std[,1:6])

# Site colors
col.site<-as.numeric(as.factor(metadata_plotspp$Site[match(labels(endo.rda.trt$CCA$u)[[1]],metadata_plotspp$Sample_ID)]))
shp.site<-c(3,8,21,22,23,24,25)[col.site]

# Plot with site as color
plot(endo.rda.trt, type='n')
points(endo.rda.trt, 'sites', pch=shp.site, cex=1.4, bg=col.site)
text(endo.rda.trt, 'bp', col='red', cex=0.9, arrow.mul=1.0)

plot(endo.rda.trt)
summary(endo.rda.trt)

ordiR2step(rda(otu.endo.otu~1, data=otu.endo.trt.std[,1:6],na.action=na.exclude), 
           scope= formula(endo.rda.trt), direction= "forward", R2scope=TRUE, pstep=1000)


##########################
# Dissimilarity analysis #
##########################

# Bray curtis dissimilarity within and among sites
# On all paired samples
meandist(vegdist(data.frame(fgh@otu_table), 'bray'),fgh@sam_data$Site) # Epi: Within groups: 0.569, between groups: 0.661    | Endo: Within groups: 0.959, between groups: 0.963
summary(meandist(vegdist(data.frame(fgh@otu_table), 'bray'),fgh@sam_data$Site))
summary(meandist(vegdist(data.frame(fgh@otu_table), 'jaccard'),fgh@sam_data$Site))

# On v3 v4
meandist(vegdist(data.frame(fgh@otu_table), 'bray'),fgh@sam_data$Site) # Epi: Within groups: 0.588, between groups: 0.654    | Endo: Within groups: 0.926, between groups: 0.952
summary(meandist(vegdist(data.frame(fgh@otu_table), 'bray'),fgh@sam_data$Site))
summary(meandist(vegdist(data.frame(fgh@otu_table), 'jaccard'),fgh@sam_data$Site))

#  Tests whether within group (but the method recommends using the anova instead)
fgh<-subset_samples(ps.pair.acesac, Life_form=='Endophyte')
fgh<-prune_taxa(taxa_sums(fgh)>0, fgh)
#diss<-mrpp(data.frame(fgh@otu_table), fgh@sam_data$Site, distance='bray')

# Partitioning beta diversity -> Cannot find GRA because total abundance is the same for each site.
beta.multi.abund(data.frame(fgh@otu_table)) # 92% bray dissimilarity for epiphytes , 98% for endophytes -> Total (includes variance within sites and among sites)
beta.multi(decostand(data.frame(fgh@otu_table),'pa')) # 92% sorensen dissimilarity for epiphytes , 97% for endophytes -> Total (include variance within sites and among sites)

# Mean within site dissimilarity

######################
## Network analyses ##
######################

#### Build networks

### Correlation matrix

# Use only 3 or 4 samples from each site to do this so every site contributes the same amount of samples
# v3 (epi) and v4 (endo) paired and 3-subsamples ps objects

# Run this using c++ implementation of sparcc algorithm (https://github.com/scwatts/fastspar): SUPER FAST vs the R version
# fastspar installable in conda env

# Overall network for epiphytes (all sites)

# Export otu table as a .tsv format
# https://forum.qiime2.org/t/exporting-otu-table-from-phyloseq-into-either-biom-or-text-format/19103/6
write_biom_tsv <- function(ps, file, sep = "; ") {
  phyloseq::t(otu_table(ps)) %>% # transpose otu matrix so asv are rows
    as.data.frame() %>%
    rownames_to_column("#OTU ID") -> phyloseq_biom
  #left_join(phyloseq::tax_table(ps) %>% 
  #            as.data.frame() %>%
  #            rownames_to_column("#OTU ID") %>% 
  #            tidyr::unite("taxonomy", !`#OTU ID`, sep = sep)) -> phyloseq_biom
  
  write_tsv(phyloseq_biom, file = file)
}

fgh3<-subset_samples(ps.pair.acesac, Life_form=='Epiphyte')
fgh3<-prune_taxa(taxa_sums(fgh3)>25, fgh3)

fgh4<-subset_samples(ps.pair.acesac, Life_form=='Endophyte')
fgh4<-prune_taxa(taxa_sums(fgh4)>25, fgh4)

write_biom_tsv(fgh3, '/home/lajoie/Garance/sparcc/otu_epi_table_fgh.tsv') # or v3
write_biom_tsv(fgh4, '/home/lajoie/Garance/sparcc/otu_endo_table_fgh.tsv') # or v4

# Within the terminal -> follow instructions from file: /home/lajoie/Garanca/sparcc/script_sparcc.txt
# Run: fastspar --otu_table ~/otu_epi_table.tsv --correlation ~/corr_epi.tsv --covariance ~/cov_epi.tsv
# Use threshold parameter (--threshold 0.2) to determine at which correlation levels pairs should be excluded
# To calculate p values of these estimates: https://github.com/scwatts/fastspar

## Epiphytes

# Import correlation matrix
cor0_1_epi<-read.csv("/home/lajoie/Garance/sparcc/median_correlation_epi_fgh.tsv", sep = "\t", row.names=1) # median_correlation_0_1_epi.tsv
cor0_1_endo<-read.csv("/home/lajoie/Garance/sparcc/median_correlation_endo_fgh.tsv", sep = "\t", row.names=1)

# Import p-vals
pval0_1_epi<-read.csv("/home/lajoie/Garance/sparcc/pvalues_epi_fgh.tsv", sep = "\t", row.names=1) # threshold_0_1_pvalues.tsv
pval0_1_endo<-read.csv("/home/lajoie/Garance/sparcc/pvalues_endo_fgh.tsv", sep = "\t", row.names=1)


# Compare thresholds at p<=0.05

pval0_1_epi[pval0_1_epi>0.05]<-NA
pval0_1_epi[pval0_1_epi<=0.05]<-1
pval0_1_epi[is.na(pval0_1_epi)==T]<-0

pval0_1_endo[pval0_1_endo>0.05]<-NA
pval0_1_endo[pval0_1_endo<=0.05]<-1
pval0_1_endo[is.na(pval0_1_endo)==T]<-0

# total number of links
pval0_1.tri_epi<-pval0_1_epi
pval0_1.tri_epi[upper.tri(pval0_1.tri_epi, diag = TRUE)]<-NA # Symmetrical matrix, need to keep only 1 half, drop the diagonal
sum(pval0_1.tri_epi, na.rm=T) # 5164 pairwise links

pval0_1.tri_endo<-pval0_1_endo
pval0_1.tri_endo[upper.tri(pval0_1.tri_endo, diag = TRUE)]<-NA # Symmetrical matrix, need to keep only 1 half, drop the diagonal
sum(pval0_1.tri_endo, na.rm=T) # 2181 pairwise links

# How many positive and negative
tot.cor_epi<-cor0_1_epi*pval0_1_epi
tot.cor.tri_epi<-tot.cor_epi
tot.cor.tri_epi[upper.tri(tot.cor.tri_epi, diag=T)]<-NA

length(tot.cor.tri_epi[tot.cor.tri_epi>0&is.na(tot.cor.tri_epi)==F]) # 3112 positive
length(tot.cor.tri_epi[tot.cor.tri_epi<0&is.na(tot.cor.tri_epi)==F]) # 2052 negative

tot.cor_endo<-cor0_1_endo*pval0_1_endo
tot.cor.tri_endo<-tot.cor_endo
tot.cor.tri_endo[upper.tri(tot.cor.tri_endo, diag=T)]<-NA

length(tot.cor.tri_endo[tot.cor.tri_endo>0&is.na(tot.cor.tri_endo)==F]) # 2180 positive
length(tot.cor.tri_endo[tot.cor.tri_endo<0&is.na(tot.cor.tri_endo)==F]) # 1 negative



# Mean links per ASV
tot.cor.pa.epi<-tot.cor_epi
tot.cor.pa.epi[tot.cor.pa.epi!=0]<-1
mean.link.epi<-rowSums(tot.cor.pa.epi)
hist(mean.link.epi)
mean(mean.link.epi) # 25 links
median(mean.link.epi) # 19

tot.cor.pa.endo<-tot.cor_endo
tot.cor.pa.endo[tot.cor.pa.endo!=0]<-1
mean.link.endo<-rowSums(tot.cor.pa.endo)
hist(mean.link.endo)
mean(mean.link.endo) # 15.5 links
median(mean.link.endo) # 16

# Make graph to characterize network
# https://kateto.net/networks-r-igraph


adj.mat<-round(as.matrix(tot.cor_epi)*100)
# adj.mat[adj.mat>-40&adj.mat<40]<-0
# Only positive links
adj.mat[adj.mat<0]<-0
adj.mat<-adj.mat[-which(rowSums(adj.mat)==0),-which(colSums(adj.mat)==0)]
dim(adj.mat)
adj.mat.epi<-adj.mat

epi.graph<-graph_from_adjacency_matrix(adj.mat.epi,
                                       mode = "undirected", weighted=TRUE,
                                       diag = FALSE)

# Color by taxonomy
col.ver<-data.frame(ps.acesac@tax_table)$class[match(names(V(epi.graph)), rownames(data.frame(ps.acesac@tax_table)))]
col.ver<-as.numeric(as.factor(col.ver))

# Color by site where the ASV is most abundant
gf<-melt(epi.tab)

gf2<-NULL
for (i in unique(gf$Var2)){
  sub<-gf[which(gf$Var2==i),]
  gf2<-rbind(gf2,sub[which(sub$value==max(sub$value)),])
}

gf3<-gf2[!duplicated(gf2[,c('Var2')]),]
gf3<-merge(gf3,metadata_plotspp[,c('Site','Sample_ID', 'Type','Plot')], by.x='Var1', by.y='Sample_ID', all.x=T)
gf3$plot<-paste(gf3$Type,gf3$Plot, sep='-')
col.ver.site<-as.numeric(as.factor(gf3$Site[match(names(V(epi.graph)), gf3$Var2)]))
#col.ver.site<-c('mediumorchid','honeydew3','yellow','tan1','deepskyblue','red','lightgreen')[col.ver.site]
shp.ver.site<-as.numeric(as.factor(gf3$plot[match(names(V(epi.graph)), gf3$Var2)]))
#shp.ver.site<-c('circle','square','triangle')[shp.ver.site]
#col.ver.site<-as.numeric(as.factor(gf3$Var1[match(names(V(epi.graph)), gf3$Var2)]))
#col.ver.site<-c("coral3","chocolate4","chartreuse4","cadetblue4","burlywood4","brown4","blue4","black","azure4","aquamarine3","deeppink4","goldenrod3","firebrick1","steelblue2","lightsalmon2","purple","plum1","khaki1","yellowgreen","ivory3","darkgreen")[col.ver.site]

# Plot
plot(epi.graph, vertex.size=4, vertex.label=shp.ver.site, vertex.color=col.ver.site, vertex.label.cex=0.7)

plot(epi.graph, vertex.size=4, vertex.label=gf3$Site[match(names(V(epi.graph)), gf3$Var2)], vertex.color=col.ver.site, vertex.label.cex=0.7)


adj.mat<-round(as.matrix(tot.cor_endo)*100)
# adj.mat[adj.mat>-40&adj.mat<40]<-0
# Only positive links
adj.mat[adj.mat<0]<-0
#adj.mat<-adj.mat[-which(rowSums(adj.mat)==0),-which(colSums(adj.mat)==0)]
dim(adj.mat)
adj.mat.endo<-adj.mat

endo.graph<-graph_from_adjacency_matrix(adj.mat.endo,
                                        mode = "undirected", weighted=TRUE,
                                        diag = FALSE)

# Color by taxonomy
col.ver<-data.frame(ps.acesac@tax_table)$class[match(names(V(endo.graph)), rownames(data.frame(ps.acesac@tax_table)))]
col.ver<-as.numeric(as.factor(col.ver))

# Color by site where the ASV is most abundant
gf<-melt(endo.tab)

gf2<-NULL
for (i in unique(gf$Var2)){
  sub<-gf[which(gf$Var2==i),]
  gf2<-rbind(gf2,sub[which(sub$value==max(sub$value)),])
}

gf3<-gf2[!duplicated(gf2[,c('Var2')]),]
gf3<-merge(gf3,metadata_plotspp[,c('Site','Sample_ID', 'Type','Plot')], by.x='Var1', by.y='Sample_ID', all.x=T)
gf3$plot<-paste(gf3$Type,gf3$Plot, sep='-')
col.ver.site<-as.numeric(as.factor(gf3$Site[match(names(V(endo.graph)), gf3$Var2)]))
#col.ver.site<-c('mediumorchid','honeydew3','yellow','tan1','deepskyblue','red','lightgreen')[col.ver.site]
shp.ver.site<-as.numeric(as.factor(gf3$plot[match(names(V(endo.graph)), gf3$Var2)]))
# By site
#col.ver.site<-as.numeric(as.factor(gf3$Site[match(names(V(endo.graph)), gf3$Var2)]))

# By sample
#col.ver.site<-as.numeric(as.factor(gf3$Var1[match(names(V(endo.graph)), gf3$Var2)]))
#col.ver.site<-c("coral3","chocolate4","chartreuse4","cadetblue4","burlywood4","brown4","blue4","black","azure4","aquamarine3","deeppink4","goldenrod3","firebrick1","steelblue2","lightsalmon2","purple","plum1","khaki1","yellowgreen","ivory3","darkgreen")[col.ver.site]

# Plot
plot(endo.graph, vertex.size=4, vertex.label=shp.ver.site, vertex.color=col.ver.site, vertex.label.cex=0.7)

plot(endo.graph, vertex.size=4, vertex.label=gf3$Site[match(names(V(endo.graph)), gf3$Var2)], vertex.color=col.ver.site, vertex.label.cex=0.7)


# Endo network super connected, despite high dissimilarity among sites


# Network characteristics

# Test the effect of grouping variables (Sample nested within site) on network structure
# https://rpubs.com/pjmurphy/338798

### Use this: https://ona-book.org/similarity.html assortativity_nominal() in igraph
# "Habitats within the plant root differ in bacterial network topology and taxonomic assortativity"
### Think about removing low-abundance taxa from the network building -> DONE

# Test if shared ASVs tend to be of certain taxonomic ID
endo.graph.net <- asNetwork(endo.graph)
plot(endo.graph.net)

endo.site<-get.vertex.attribute(endo.graph2, 'Site')
endo.sample<-get.vertex.attribute(endo.graph2, 'Sample_ID')

# Plot
#https://cran.r-project.org/web/packages/bipartite/vignettes/Intro2bipartite.pdf

## Endophyte network

# Modularity -> Endophyte highly modular, driven by uniqueness of samples (ASVs occur in only 1 sample, and thus they are 'linked')
modules.endo=computeModules(adj.mat.endo)
plotModuleWeb(modules.endo)
mod.id<-listModuleInformation(modules.endo)
mm<-fgh4@otu_table[,which(colnames(fgh4@otu_table)%in%unlist(mod.id[[2]][[2]]))]

# Calculate modularity
men<-DIRT_LPA_wb_plus(modules.endo@modules) # 0.45

# Calculate nestedness
nest.endo<-nest.smdm(adj.mat.endo, weighted=T)
const<-module2constraints(modules.endo)
nest.endo.mod<-nest.smdm(adj.mat.endo, constraint=const, weighted=T)

# Calculate H2
spec.endo<-H2fun(adj.mat.endo) # 0.48


# Components
# Calculate the maximal (weakly or strongly) connected components of a graph
compon.endo<-components(graph.endo)
groups(compon.endo)
# Corresponds ~ to first layer of modularity -> Could be focusing epiphyte analyses on the biggest component to decrease complexity

# Calculate connectedness
conn.endo<-edge_density(endo.graph) # 1

# Transitivity
# Transitivity measures the probability that the adjacent vertices of a vertex are connected. This is sometimes also called the clustering coefficient.
#trans.endo<-transitivity(endo.graph) # 0.68

# Testing network properties by bootstrapping -> check R script 'networkbootstraps.R'
load("/home/lajoie/Garance/boot_endo_objects.RData")

# Connectance
hist(conn.endo.boot)
abline(v=conn.endo, col='red')
length(conn.endo.boot[conn.endo.boot>conn.endo])/1000 # 0.17

# Modularity
hist(modul.endo.boot)
abline(v=men$modularity, col='red')
length(modul.endo.boot[modul.endo.boot<men$modularity])/1000 # 0.061

# Nestedness
hist(nested.endo.boot)
abline(v=nest.endo$WNODFmatrix, col='red')
length(nested.endo.boot[nested.endo.boot>nest.endo$WNODFmatrix])/1000 # 0.126

# Specialization
hist(h2.endo.boot)
abline(v=spec.endo['H2'], col='red')
length(h2.endo.boot[h2.endo.boot<spec.endo['H2']])/1000 # 0.035

# Prop positive link
prop.endo.posi<-pos.link.endo.boot/(pos.link.endo.boot+neg.link.endo.boot)
hist(prop.endo.posi)
abline(v=2180/2181, col='red')
length(prop.endo.posi[prop.endo.posi>(2180/2181)])/1000 # 0.015

# Assortativity -> Needs to be calculated based on some characteristic of ASVs (e.g. site where most abundant)
V(endo.graph)$Site<-gf3$Site[match(names(V(endo.graph)),gf3$Var2)]
V(endo.graph)$Plot<-gf3$Var1[match(names(V(endo.graph)),gf3$Var2)]

V(endo.graph)$class<-tax_table(ps.acesac)[,3][match(names(V(endo.graph)),rownames(tax_table(ps.acesac)[,3]))]
V(endo.graph)$class[is.na(V(endo.graph)$class)==T]<-'Unknown'

an.site<-igraph::assortativity_nominal(
  endo.graph, 
  as.integer(as.factor(V(endo.graph)$Site))
) # 0.68

# Testing assortativity significance by random attribution of properties

boots<-NULL

for (i in 1:999){
  
  V(endo.graph)$site.rand<-sample(V(endo.graph)$Site)
  
  bn<-igraph::assortativity_nominal(
    endo.graph, 
    as.numeric(as.factor(V(endo.graph)$site.rand))
  ) 
  boots<-c(boots,bn)
  
}

hist(boots)

# p-val
length(boots[boots>an.site])/length(boots)


an.plot<-igraph::assortativity_nominal(
  endo.graph, 
  as.integer(as.factor(V(endo.graph)$Plot))
) # 0.68

boots<-NULL

for (i in 1:999){
  
  V(endo.graph)$plot.rand<-sample(V(endo.graph)$Plot)
  
  bn<-igraph::assortativity_nominal(
    endo.graph, 
    as.numeric(as.factor(V(endo.graph)$plot.rand))
  ) 
  boots<-c(boots,bn)
  
}

hist(boots)

# p-val
length(boots[boots>an.plot])/length(boots)



an.class<-igraph::assortativity_nominal(
  endo.graph, 
  as.numeric(as.factor(V(endo.graph)$class))
) # 0.68

# Testing assortativity significance by random attribution of properties

boots<-NULL

for (i in 1:999){
  
  V(endo.graph)$class.rand<-sample(V(endo.graph)$class)
  
  bn<-igraph::assortativity_nominal(
    endo.graph, 
    as.numeric(as.factor(V(endo.graph)$class.rand))
  ) 
  boots<-c(boots,bn)
  
}

hist(boots)

# p-val
length(boots[boots>an.class])/length(boots)


# Epiphyte

modules.epi=computeModules(adj.mat.epi, forceLPA=TRUE) # Need to use force argument otherwise RStudio crashes
#plotModuleWeb(modules.epi)
mod.id2<-listModuleInformation(modules.epi)
mm<-v3@otu_table[,which(colnames(v3@otu_table)%in%unlist(mod.id2[[2]][[4]]))]

# Calculate modularity
mep<-DIRT_LPA_wb_plus(modules.epi@modules) # 0.40
# Calculate nestedness
nest.epi<-nest.smdm(adj.mat.epi, weighted=T)
const.epi<-module2constraints(modules.epi)
nest.epi.mod<-nest.smdm(adj.mat.epi, constraint=const.epi, weighted=T)

# Calculate H2
spec.epi<-H2fun(adj.mat.epi) # 0.45

# Calculate connectance
conn.epi<-edge_density(epi.graph) # 1

# Transitivity
# Transitivity measures the probability that the adjacent vertices of a vertex are connected. This is sometimes also called the clustering coefficient.
#trans.epi<-transitivity(epi.graph) # 0.41

# Testing network properties by bootstrapping -> check R script 'networkbootstraps.R'
#load("/home/lajoie/Garance/boot_endo_objects.RData")
load("/home/lajoie/Garance/boot_epi_objects.RData")

# Connectance
hist(conn.epi.boot)
abline(v=conn.epi, col='red') # Out of graph, higher
length(conn.epi.boot[conn.epi.boot>conn.epi])/1000 # 0

# Modularity
hist(modul.epi.boot)
abline(v=mep$modularity, col='red') # Out of graph, lower
length(modul.epi.boot[modul.epi.boot<mep$modularity])/1000 # 0

# Nestedness
hist(nested.epi.boot)
abline(v=nest.epi$WNODFmatrix, col='red') # Out of graph, higher
length(nested.epi.boot[nested.epi.boot>nest.epi$WNODFmatrix])/1000 # 0

# Specialization
hist(h2.epi.boot)
abline(v=spec.epi['H2'], col='red') # Out of graph, lower
length(h2.epi.boot[h2.epi.boot<spec.epi['H2']])/1000 # 0

# Prop positive link
prop.epi.posi<-pos.link.epi.boot/(pos.link.epi.boot+neg.link.epi.boot)
hist(prop.epi.posi)
abline(v=3112/5164, col='red') # Out of graph, higher
length(prop.epi.posi[prop.epi.posi>(3112/5164)])/1000 # 0.015

# Assortativity -> Needs to be calculated based on some characteristic of ASVs (e.g. site where most abundant)
# Is a measure of the extent to which vertices with the same properties connect to each other.

V(epi.graph)$Site<-gf3$Site[match(names(V(epi.graph)),gf3$Var2)]
V(epi.graph)$Plot<-gf3$Var1[match(names(V(epi.graph)),gf3$Var2)]

V(epi.graph)$class<-tax_table(ps.acesac)[,3][match(names(V(epi.graph)),rownames(tax_table(ps.acesac)[,3]))]
V(epi.graph)$class[is.na(V(epi.graph)$class)==T]<-'Unknown'


an.site<-igraph::assortativity_nominal(
  epi.graph, 
  as.integer(as.factor(V(epi.graph)$Site))
) # 0.31

# Testing assortativity significance by random attribution of properties

boots<-NULL

for (i in 1:999){
  
  V(epi.graph)$site.rand<-sample(V(epi.graph)$Site)
  
  bn<-igraph::assortativity_nominal(
    epi.graph, 
    as.numeric(as.factor(V(epi.graph)$site.rand))
  ) 
  boots<-c(boots,bn)
  
}

hist(boots)

# p-val
length(boots[boots>an.site])/(length(boots)+1)


##
an.plot<-igraph::assortativity_nominal(
  epi.graph, 
  as.integer(as.factor(V(epi.graph)$Plot))
) # 0.24

boots<-NULL

for (i in 1:999){
  
  V(epi.graph)$plot.rand<-sample(V(epi.graph)$Plot)
  
  bn<-igraph::assortativity_nominal(
    epi.graph, 
    as.numeric(as.factor(V(epi.graph)$plot.rand))
  ) 
  boots<-c(boots,bn)
  
}

hist(boots)

# p-val
length(boots[boots>an.plot])/length(boots)


an.class<-igraph::assortativity_nominal(
  epi.graph, 
  as.numeric(as.factor(V(epi.graph)$class))
) # 0.33-

# Testing assortativity significance by random attribution of properties

boots<-NULL

for (i in 1:999){

V(epi.graph)$class.rand<-sample(V(epi.graph)$class)

bn<-igraph::assortativity_nominal(
  epi.graph, 
  as.numeric(as.factor(V(epi.graph)$class.rand))
) 
boots<-c(boots,bn)

}

hist(boots)

# p-val
length(boots[boots>an.class])/length(boots)

