#######################################
### C. sakazakii Pangenome Analysis ###
#######################################

# Mairui Gao
# mrgao@umd.edu
# last edited:07/02/2025

######################################

## set directory
setwd("~/Cronobacter")

library(ape)
library(caret)
library(DMwR)
library(dplyr)
library(dunn.test)
library(ggalluvial)
library(ggplot2)
library(minpack.lm)
library(patchwork)
library(phangorn)
library(pROC)
library(reshape2)
library(tidyr)
library(vegan)



### pangenome
pangenome_mat <- read.table("Pan_analysis/gene_presence_absence.Rtab", h=T, row.names=1)
meta <- read.csv("metadata/metadata.csv", h=T, row.names=1)

# core and accessory
genes_by_genome <- apply(pangenome_mat, 1, sum)
gene_freq1 <- data.frame(genome = c(1:dim(pangenome_mat)[2]),
                       count = c(1:dim(pangenome_mat)[2]),
                       group = c(
                         rep("Accessory", floor(dim(pangenome_mat)[2]*0.95)),
                         rep("Core", ceiling(dim(pangenome_mat)[2]*0.05))))

# core, shell and cloud
# assign the first 15% to cloud, the next 80% to shell, and the last 5% to core
gene_freq2 <- data.frame(genome = c(1:dim(pangenome_mat)[2]),
                       count = c(1:dim(pangenome_mat)[2]),
                       group = c(rep("Cloud", floor(dim(pangenome_mat)[2]*0.15)),
                                 rep("Shell", ceiling(dim(pangenome_mat)[2]*0.80)),
                                 rep("Core", round(dim(pangenome_mat)[2]*0.05))))


gene_freq2$group <- factor(gene_freq2$group, levels = c("Core", "Shell","Cloud"))

# calculate the number of genes are present in i genomes and stores in gene_freq
for (i in 1:dim(pangenome_mat)[2]) {
  gene_freq2[i,2] = length(which(genes_by_genome == i))
}

ggplot(gene_freq2, aes(x = genome, y = count, fill = group)) +
  geom_bar(stat = "identity", width = 4#, fill = "red"
  ) +
  xlab("Number of genomes") +
  ylab("Gene clusters") +
  theme_classic() +
  #scale_y_log10() +
  xlim(-10,800) +
  scale_fill_manual(values=c("#FB9A99","gray","#80B1D3")) +
  theme(axis.text = element_text(size=16, color = "black"), 
        axis.title = element_text(size=18, color = "black")) +
  theme(legend.position="none")

ggsave('gene_freq2.pdf', width=5, height=4)

pie.df <- data.frame(count = tapply(gene_freq2$count, gene_freq2$group, sum),
                     percent = 100*tapply(gene_freq2$count, gene_freq2$group, sum)/sum(gene_freq2$count))
pie.df$percent = round(pie.df$percent, 1)
pie.df$group<-rownames(pie.df)

ggplot(pie.df, aes(x = "", y = percent, fill = group)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Core" = "#FB9A99", "Shell" = "gray", "Cloud" = "#80B1D3")) +
  geom_text(aes(label = paste0(round(percent, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 6) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16))

ggsave('gene_distribution_pie.pdf', width=5, height=4)

# identify gene type 
total_genomes <- dim(pangenome_mat)[2]
core_threshold <- ceiling(total_genomes * 0.95)

gene_type <- data.frame(genes_by_genome,
                        category = rep("Core", length(genes_by_genome)))
gene_type$category <- as.character(gene_type$category)
gene_type$category[gene_type$genes_by_genome <= core_threshold] <- "Accessory"

# Determine fraction count per genome                              
genes_per_genome <- matrix(nrow = dim(pangenome_mat)[2], ncol = 2)

for (i in 1:dim(pangenome_mat)[2]) {
  genes_per_genome[i,] = rbind(tapply(pangenome_mat[,i], gene_type$category, sum))
}

rownames(genes_per_genome) <- colnames(pangenome_mat)
colnames(genes_per_genome) <- names(tapply(pangenome_mat[,1], gene_type$category, sum))


# accessory and core genome
core<-subset.data.frame(gene_type, category=="Core")
core$Variable<-row.names(core)
core<-core[,-1]
core<-core[,-1]
core<-as.data.frame(core)
colnames(core)[1]<-c('Variable')
pangenome_mat$Variable<-row.names(pangenome_mat)

dfcore<- pangenome_mat %>% 
  semi_join(core, by = "Variable")
dfcore<-dfcore[,-749]
pangenome_mat<-pangenome_mat[,-749]

# check for colnames and subset
colnames(dfcore) == colnames(pangenome_mat)
dfaccessory <- anti_join(pangenome_mat, dfcore, by = colnames(dfcore))

write.csv(dfcore,'dfcore.csv')
write.csv(dfaccessory,'dfaccessory.csv')

### permutation
## skip to 'load result' if already completed permutation

## gene count function counts how many genes are first found in each genome (column)
new_gene_calc <- function (x) {
  gene = c(1:dim(x)[1])
  for (i in 1:dim(x)[1]) {
    gene[i] = which(x[i,] > 0)[1] # Find the first genome with the gene
  }
  gene_count = c(1:dim(x)[2]) 
  for (z in 1:dim(x)[2]) {
    gene_count[z] = length(which(gene == z)) # Count genes first seen in genome z
  }
  print(gene_count)
}

## 100x permutations
perm_pan_bac <- as.data.frame(matrix(nrow = 100, ncol = dim(pangenome_mat)[2]))
for (n in 1:100) {
  perm_pan_bac[n,] = new_gene_calc(pangenome_mat[,sample(dim(pangenome_mat)[2])])
  print(n)
} 
write.table(perm_pan_bac, 'permutation_new_genes_per_genome.txt', sep="\t", col.names=NA, quote = FALSE)

## load result (instead of re-running)
perm_pan <- read.table('permutation_new_genes_per_genome.txt', h=T, row.names=1)

## model and plot: set data
newgenes = perm_pan
colnames(newgenes) = c(1:dim(newgenes)[2])
newgenes_m = melt(newgenes)
head(newgenes_m)
newgenes_stat = data.frame(genomes = c(1:dim(pangenome_mat)[2]),
                           ng_mean = tapply(newgenes_m$value, newgenes_m$variable, mean),
                           ng_median = tapply(newgenes_m$value, newgenes_m$variable, median),
                           ng_sd = tapply(newgenes_m$value, newgenes_m$variable, sd),
                           ng_se = tapply(newgenes_m$value, newgenes_m$variable, sd)/sqrt(100))

tapply(newgenes_m$value, newgenes_m$variable, mean)

# power-law fit 
ng_power = nlsLM(log(ng_mean+1)~a*genomes^-b, 
                 data=newgenes_stat, start = list(a=0.1,b=0.1))
summary(ng_power)  
ng_power_predict = data.frame(var = 1:dim(pangenome_mat)[2],
                              val = predict(ng_power))
# N100 value
exp(ng_power_predict$val[100])

## plot model fit to averages, displaying all permuted data points
ggplot(data = newgenes_m, aes(x = as.numeric(variable), y = value+1)) +
  geom_point(alpha = 0.1, shape = 1, col = "gray", size = 0.7) +
  #geom_point(data = newgenes_stat, aes(x = genomes, y = ng_mean+1), shape = 1, size = 1.5, alpha=0.8, col = "black") +
  geom_line(data = newgenes_stat, aes(x = genomes, y = ng_mean+1), col = "blue", size=0.8) +
  #geom_line(data = ng_power_predict, aes(x = var, y = exp(val)), lty = 1, col = "red", size=0.8) +
  theme_classic() +
  xlab("Number of genomes") +
  ylab("New gene clusters") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500, 600,700), limits = c(0,760)) +
  scale_y_log10(breaks=c(1,4,10,40,100,400,1000,4000), limits = c(1,8000)) +
  theme(axis.text = element_text(size=16, color = "black"), 
        axis.title = element_text(size=18, color = "black"))

ggsave('permutation.pdf', width=5, height=4)


# Mantel test for correlation with accessory gene diversity
# read tree
core_tree_bac <- read.tree("Pan_analysis/mlst.nwk")
pangenome_bac <- read.table("Pan_analysis/gene_presence_absence.Rtab", h=T, row.names=1)

# prep core gene distance matrix
bac_core_distance <- cophenetic.phylo(midpoint(core_tree_bac))
bac_core_distance <- bac_core_distance[order(rownames(bac_core_distance)), order(colnames(bac_core_distance))]

# prep accessory gene distance matrix
bac_accessory_distance <- vegdist(t(pangenome_bac[which(apply(pangenome_bac, 1, sum)<0.95*dim(pangenome_bac)[2]),]), 
                                  'jaccard', binary = T)
bac_accessory_distance <- as.matrix(bac_accessory_distance)
bac_accessory_distance <- bac_accessory_distance[order(rownames(bac_accessory_distance)), order(colnames(bac_accessory_distance))]

# check for same order
colnames(bac_core_distance) == colnames(bac_accessory_distance)
rownames(bac_core_distance) == rownames(bac_accessory_distance)

# run mantel test
mantel(bac_core_distance, bac_accessory_distance)


## metadata sankey diagram
sankey <- read.csv("metadata/meta_sankey.csv")
ggplot(data = sankey, 
       aes(axis1 = Sample_type_specific, axis2 = Sample_type, axis3 = Continent, y = Total)) +
       geom_alluvium(aes(fill = Sample_type_specific), width = 0.1, alpha = 0.8) +  # Flow color
       geom_stratum(aes(fill = after_stat(stratum)), width = 0.1, alpha = 0.8) +  # Node color
       geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 0, size = 4) +
       scale_x_discrete(limits = c("Sample_type_specific", "Sample_type", "Continent")) +
       scale_fill_manual(values = c(
                         "Clin_Blood+CSF" = "#BC80BD","Clin_Other" = "#CAB2D6","Clin_Stool" = "#810F7C",
                         "Env_Built" = "#66C2A5","Env_FoodFacility" = "#B3DE69","Env_Natural" = "#238443",
                         "Food_Dry" = "#FB6A4A","Food_Other" = "#E41A1C","Food_Powdered" = "#FB9A99",
                         "Clinical" = "#CAB2D6","Environment" = "#7FC97F","Food" = "#FBB4AE",  
                         "Asia" = "#E5D8BD","Europe" = "#DFC27D","North_America" = "#BF812D")) +
      theme_bw() +
      theme(
            legend.position = 'right',
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14),
            axis.title = element_blank(),
            axis.text.x = element_text(face = 'bold', colour = 'black', size = 12),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
  )

ggsave('meta_sankey.pdf', width=12, height=6)



### beta diversity on accessory genome
meta <- read.csv("metadata/metadata.csv", h=T)
meta <- subset.data.frame(meta, Continent!="South_America")

# dfa was generated previously
dfa <- read.csv('dfaccessory.csv', row.names=1, h=T)
dfa <- dfa[,meta$Bin]

# for subsets
# meta.a<-subset.data.frame(meta, Continent =="Asia")
# df1<-dfa[,meta.a$Bin]

# subsets df1 may have rows which contain absent genes (all 0)
# gene_absence <- rowSums(df1 != 0) == 0
# df1 <- df1[!gene_absence,]

# continent 
# compare columns
beta <- vegdist(t(dfa), method = 'jaccard', binary = TRUE)
# beta dispersion, similar to assess variance in ANOVA
beta.mat <- as.matrix(beta)
pcoa <- cmdscale(beta, k = 4, eig = TRUE)
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') 
ord$Bin = rownames(ord)
beta_df<-merge(ord,meta, by="Bin")

# beta dispersion and PERMDISP test
#dispersion <- betadisper(beta, beta_df$Sample_type)
dispersion <- betadisper(beta, beta_df$Continent)
print(dispersion)

permutest(dispersion, permutations = 999)
plot(dispersion)


### PCoA with side boxplot
# get percentage of pcoa1 and pcoa2 
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# make Continent/Source as group
colnames(beta_df)[6]<-c('group')
beta_df$group<-factor(beta_df$group)

mycol1<-c("Asia" = "#E5D8BD", 
          "Europe" = "#DFC27D", 
          "North_America" = "#BF812D")

mycol2<-c("Clinical" = "#CAB2D6", 
          "Environment" = "#7FC97F", 
          "Food" = "#FBB4AE")

mycol3<-c("Food_Dry" = "#FB6A4A", 
          "Food_Other" = "#E41A1C", 
          "Food_Powdered" = "#FB9A99")

p <- ggplot(beta_df, aes(x = pcoa1, y = pcoa2, color = group)) +
  geom_point(size = 3, shape=17) +  # shape=16, 18
  scale_color_manual(values = mycol1) +  # mycol2 # mycol3
  geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "solid") +
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "solid") +
  theme_classic() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  
  )

df1 <- beta_df %>%
  mutate(group = recode(group,
                        "Asia" = "A",
                        "Europe" = "E",
                        "North_America" = "N"))

mycolc<-c("A" = "#E5D8BD", 
          "E" = "#DFC27D", 
          "N" = "#BF812D")

df2 <- beta_df %>%
  mutate(group = recode(group,
                        "Clinical" = "C",
                        "Environment" = "E",
                        "Food" = "F"))

mycols<-c("C" = "#CAB2D6", 
          "E" = "#7FC97F", 
          "F" = "#FBB4AE")

df3 <- beta_df %>%
  mutate(group = recode(group,
                        "Food_Dry" = "FD",
                        "Food_Other" = "FO",
                        "Food_Powdered" = "FP"))

mycolf<-c("FD" = "#FB6A4A", 
          "FO" = "#E41A1C", 
          "FP" = "#FB9A99")

pb <-ggplot(df1, aes(x = pcoa1, y=group, fill=group)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.4) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = mycolc) +  # mycols # mycolf
  xlab("PCo1(%)") +
  theme_classic() +  
  theme(
    legend.position = "none",  
    axis.text.y = element_text(size = 16, angle = 0),  
    axis.text.x = element_text(size = 16, angle = 0),  
    axis.title.y = element_blank(),  
    axis.title.x = element_text(size = 20)
  )

pl<-ggplot(df1, aes(y = pcoa2, x=group, fill=group)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.4) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = mycolc) + 
  ylab("PCo2(%)") +
  theme_classic() +  
  theme(
    legend.position = "none",  
    axis.text.y = element_text(size = 16, angle = 0),  
    axis.text.x = element_text(size = 16, angle = 0),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 20)
  )


p0 <- ggplot() + theme(panel.background = element_blank())

pl+p+p0+pb+plot_layout(ncol = 2,nrow = 2,heights = c(4,1),widths = c(1,5))& 
  theme(plot.margin = margin(0, 0, 0, 0)) # remove margin

ggsave('_pcoa.pdf', width=6, height=6)


## centroids plot
cen1 <- beta_df %>%
  group_by(group) %>%
  summarize(
    se_pcoa1 = sd(pcoa1) / sqrt(n()),  
    se_pcoa2 = sd(pcoa2) / sqrt(n()))

cen2 <- beta_df %>%
  group_by(group) %>%
  summarize(
    pcoa1 = mean(pcoa1),
    pcoa2 = mean(pcoa2))

cen<-cbind(cen1, cen2)
cen<-cen[,-4]
colnames(cen)[1]<-c('var')

ggplot(cen, aes(x = pcoa1, y = pcoa2, color = var)) +
  geom_point(size = 5) +
  # width in geom_errorbar and height in geom_errorbarh are defined in coordinates, need manual adjust
  # Vertical error bars 
  geom_errorbar(aes(ymin = pcoa2 - se_pcoa2, ymax = pcoa2 + se_pcoa2), width = 0.004, size = 1.5) + 
  # Horizontal error bars
  geom_errorbarh(aes(xmin = pcoa1 - se_pcoa1, xmax = pcoa1 + se_pcoa1), height = 0.008, size = 1.5) +  
  scale_color_manual(values = mycol1) +
  xlab("PCo1(%)") +
  ylab("PCo2(%)") +
  theme_classic() + 
  theme(
    axis.title = element_text(size=14),
    axis.text = element_text(size=12), 
    legend.position = "top",
    legend.title = element_blank(),
    text = element_text(size = 12)
  )

ggsave('source_centroid.pdf', width=4, height=4.2)


# pairwise PERMANOVA
# dfa was generated by previous steps
df <- read.csv('dfaccessory.csv',header = T,row.names = 1)
groups <- read.csv('metadata/metadata.csv',header = T)
groups <- subset.data.frame(groups, Continent!="South_America")
colnames(groups)[1]<-c('sample')
colnames(groups)[2]<-c('group')

# for subsets
# groups<-subset.data.frame(groups, Sample_type == "Food")
# groups<-subset.data.frame(groups, Sample_type == "Environment")
# groups<-subset.data.frame(groups, Sample_type == "Clinical")

# after subsetting the df, be careful there could be columns containing all 0, remove them
# all_zero_columns <- colSums(df != 0) == 0
# df <- df[, !all_zero_columns]

groups<-groups %>%mutate(group=factor(group,levels = unique(group)))
mycompare<-as.data.frame(combn(levels(groups$group),2))
myresults<-NULL

# make sure df has rows as genes, and columns as samples
for(i in 1:ncol(mycompare)){
  #i=1
  groups1<-groups %>% filter(group %in% mycompare[,i] ) %>% mutate(group=factor(group,levels = unique(group)))
  df1<-df[,groups1$sample]
  tdf<-as.data.frame(t(df1))
  tdf <- as.data.frame(lapply(tdf, function(x) as.numeric(as.character(x))))
  
  myadonis<-adonis2(tdf~groups1$group,permutations = 999, method='jaccard')
  ADONIS_F<-myadonis$F[1]
  ADONIS_R2<-myadonis$R2[1]
  ADONIS_P<-myadonis$`Pr(>F)`[1]
  
  myresult1<-as.data.frame(cbind(myset=paste(mycompare[,i],collapse = ' vs '),
                                 ADONIS_F,  ADONIS_R2,ADONIS_P))
  myresults<-rbind(myresult1,myresults)
}

write.csv(myresults,'df_adonis.csv')


# kruskal wallis on PCo1 and PCo2
kruskal.test(pcoa1~Sample_type, beta_df)
dunn.test(x = beta_df$pcoa1, g = beta_df$Sample_type, method = "bh")  


##  genes per genome
df<-read.table("Pan_analysis/gene_presence_absence.Rtab", h=T, row.names=1)
df1<-as.data.frame(apply(df, 2, sum))
colnames(df1)[1]<-c('sum')
meta <- read.csv("metadata.csv", h=T)
meta <- subset.data.frame(meta, Continent!="South_America")
df1$Bin<-row.names(df1)
df2<-merge(df1, meta, by="Bin")
# identify group factor for comparison
colnames(df2)[3]<-c('group')

se <- df2 %>%
  group_by(group) %>%
  summarize(
    se = sd(sum) / sqrt(n()),  
  )

mean <- df2 %>%
  group_by(group) %>%
  summarize(
    mean = mean(sum),
  )

df3<-cbind(se, mean)
df3<-df3[,-3]
df3$Lower <- df3$mean - df3$se
df3$Upper <- df3$mean + df3$se


ggplot(df3, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +  
  ylab("Genes per genome") +
  xlab("#General source") +
  # facet_wrap(~ Sample_type, scales = "free_x", ncol = 3)+
  theme_bw() +
  scale_fill_manual(values = mycol1)+
  coord_cartesian(ylim = c(3500, 4200)) + 
  theme(
    axis.title = element_text(size = 14),        
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(),
    legend.position = "right"
  )

ggsave('genes_per_genome.pdf', width=5.5, height=4)
ggsave('genes_per_genome_continent.pdf', width=7, height=4)   

kruskal.test(sum~group, df2)
dunn.test(x = df2$sum, g = df2$group, method = "bh")

# subset df for specific source comparison
df_clin<-subset.data.frame(df2, Sample_type == "Clinical")
df_env<-subset.data.frame(df2, Sample_type == "Environment")
df_food<-subset.data.frame(df2, Sample_type == "Food")


#COG
df <- read.table("Pan_analysis/gene_presence_absence.Rtab", h=T, row.names=1)
# dfaccessory.csv was generated by previous steps
# df <- read.csv('dfaccessory.csv', row.names=1, h=T)
row.names(df) <- gsub("~~~", "_", row.names(df))
cog<-read.csv('Features/cog_annotation.csv', h=T)

# split multi-categories
cog1<- cog %>%
  separate_rows(COG, sep = "") %>%
  filter(COG != "")

cog1<-cog1[,-2]
df$query<-row.names(df)
df1<-merge(df, cog1, by="query")


# make df1 numeric except for COG
num <- sapply(df1, is.numeric)
num["COG"] <- TRUE

category_sums <- aggregate(. ~ COG, df1[, num], sum)

df2<-as.data.frame(t(category_sums))
colnames(df2)<-df2[1,]
df2<-df2[-1,]
meta <- read.csv("metadata/metadata.csv", h=T)
df2$Bin<-row.names(df2)
df3<-merge(df2, meta, by="Bin")
write.csv(df3, 'cog_dfa_meta.csv')


df1<-read.csv("cog_dfa_meta.csv")
#  continent
### Fisher's Exact Test
# make df4 numeric except for COG/Sample_type
num <- sapply(df4, is.numeric)
num["Continent"] <- TRUE
# num["Sample_type"] <- TRUE
# num["Sample_type_specific"] <- TRUE

category_sums <- aggregate(. ~ Continent, df4[, num], sum)
# category_sums <- aggregate(. ~ Sample_type, df1[, num], sum)
# category_sums <- aggregate(. ~ Sample_type_specific, df1[, num], sum)

row.names(category_sums)<-category_sums[,1]
category_sums<-category_sums[,-1]
# remove South_America if necessary
category_sums<-category_sums[-4,]

# prepare contingency tables for fisher test
contingency_tables <- list()

for (cog in colnames(category_sums)) {
  contingency_table <- matrix(NA, nrow = 2, ncol = 3)
  rownames(contingency_table) <- c(cog, paste0("Non-", cog))
  colnames(contingency_table) <- rownames(category_sums)
  for (source in rownames(category_sums)) {
    contingency_table[1, source] <- category_sums[source, cog]  
    contingency_table[2, source] <- sum(category_sums[source, ]) - category_sums[source, cog]  
  }
  contingency_tables[[cog]] <- contingency_table
}

## Fisher’s test
results <- data.frame(
  COG = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each contingency table and perform Fisher's test with simulation
for (cog in names(contingency_tables)) {
  fisher_result <- fisher.test(contingency_tables[[cog]], simulate.p.value = TRUE, B = 1e5)
  results <- rbind(results, data.frame(
    COG = cog,
    P_value = fisher_result$p.value
  ))
}

results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

write.csv(results, 'results.csv')


### enriched gene
#dfa<-read.csv('dfaccessory.csv', row.names=1, h=T)
df<-read.table("Pan_analysis/gene_presence_absence.Rtab", h=T, row.names=1)
row.names(df) <- gsub("~~~", "_", row.names(df))
cog<-read.csv('Features/cog_annotation.csv')

# split multi-categories
cog1<- cog %>%
  separate_rows(COG, sep = "") %>%
  filter(COG != "")

### run every enriched COG category
# change COG L to others such as E, G, I...
cogl<-subset.data.frame(cog1, COG=="L")
cogl<-cogl[,-2]
cogl<-as.data.frame(t(cogl))
colnames(cogl)<-NULL
colnames(cogl)<-cogl[1,]
cogl<-cogl[-1,]
cogl<-cogl[-1,]

#dfa colnames as gene name
dfa<-as.data.frame(t(df))

# subset dfa with L genes
dfa_cog <- intersect(colnames(dfa), colnames(cogl))
dfcogl <- dfa[, dfa_cog, drop = FALSE]

meta <- read.csv("metadata/metadata.csv", h=T)
dfcogl$Bin<-rownames(dfcogl)
df1<-merge(dfcogl, meta, by ="Bin")

# revise df1 about sample/continent before next step
row.names(df1)<-df1[,1]
df1<-df1[,-1]
# ...

# make df1 numeric except for Continent/Sample_type/Sample_type_specific
num <- sapply(df1, is.numeric)
num["Continent"] <- TRUE
# num["Sample_type"] <- TRUE
# num["Sample_type_specific"] <- TRUE

category_sums <- aggregate(. ~ Continent, df1[, num], sum)
#category_sums <- aggregate(. ~ Sample_type, df1[, num], sum)
#category_sums <- aggregate(. ~ Sample_type_specific, df1[, num], sum)

row.names(category_sums)<-category_sums[,1]
category_sums<-category_sums[,-1]
# remove South_America
category_sums<-category_sums[-4,]
write.csv(category_sums, 'cogL_sums.csv')

# remove all 0 columns
all_zero_columns <- colSums(category_sums != 0) == 0
category_sums <- category_sums[, !all_zero_columns]


# prepare contingency tables for fisher test
contingency_tables <- list()

for (cog in colnames(category_sums)) {
  contingency_table <- matrix(NA, nrow = 2, ncol = 3)
  rownames(contingency_table) <- c(cog, paste0("Non-", cog))
  colnames(contingency_table) <- rownames(category_sums)
  for (source in rownames(category_sums)) {
    contingency_table[1, source] <- category_sums[source, cog]  
    contingency_table[2, source] <- sum(category_sums[source, ]) - category_sums[source, cog]  
  }
  contingency_tables[[cog]] <- contingency_table
}

## Fisher’s test
results <- data.frame(
  COG = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each contingency table and perform Fisher's test with simulation
for (cog in names(contingency_tables)) {
  fisher_result <- fisher.test(contingency_tables[[cog]], simulate.p.value = TRUE, B = 1e5)
  results <- rbind(results, data.frame(
    COG = cog,
    P_value = fisher_result$p.value
  ))
}

results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

write.csv(results, 'results.csv')


# cog relative abundance
# df0 <- read.csv('cog_pan.csv', row.names=1, h=T)
# cog_dfa was generated by combing dfa and cog annotation
df0 <- read.csv('cog_dfa.csv', row.names=1, h=T)
df1 <- sweep(df0,1,rowSums(df0), '/')
df <- t(df1)
df <- as.data.frame(df)
meta <- read.csv("metadata/metadata.csv", h=T)

# cog average value 
outord<-meta$Sample_type[!duplicated(meta$Sample_type)]
# meta<-subset.data.frame(meta, Continent!="South_America")
# outord<-meta$Continent[!duplicated(meta$Continent)]
# meta<-subset.data.frame(meta, Sample_type=="Food")
# outord<-meta$Sample_type_specific[!duplicated(meta$Sample_type_specific)]
df<-df[,meta$Bin]
colnames(df)<-meta$Sample_type
# colnames(df)<-meta$Continent
# colnames(df)<-meta$Sample_type_specific
df1<-as.data.frame(t(apply(df, 1, function(x) tapply(x, colnames(df), mean))))
df1<-df1[,outord,drop=F]
df2<-cbind(c(rownames(df1)),df1)
colnames(df2)[1]<-'cog'

df3<-melt(df2,id.vars = 'cog')

colset<-read.csv('Features/cog_color.csv')
colnames(colset)[1:2]<-c('variable','Color')
mycol<- colset$Color
colnames(df3)[1]<-'variable'
colnames(df3)[2]<-'Sample_type' #'Continent' 

custom_labels <- c(
  "A-RNA processing and modification",
  "B-Chromatin structue and dynamics",
  "C-Energy production and conversion",
  "D-Cell cycle control, cell division, chromosome partition",
  "E-Amino acid transport and metabolism",
  "F-Nucleotide transport and metabolism",
  "G-Carbohydrate transport and metabolism",
  "H-Coenzyme transport and metabolism",
  "I-Lipid transport and metabolism",
  "J-Translation, ribosomal structure and biogenesis",
  "K-Transcription",
  "L-Replication, recombination and repair",
  "M-Cell wall/membrane/envelop biogenesis",
  "N-Cell motility",
  "O-Posttranslational modification, protein turnover, chaperones",
  "P-Ion transport and metabolism",
  "Q-Secondary metabolites biosynthesis, transport and catabolism",
  "S-Function unknown",
  "T-Signal transduction mechanisms",
  "U-Intracellular trafficking, secretion, and vesicular transport",
  "V-Defense mechanisms"
)

ggplot(df3, aes(x=Sample_type,stratum=variable, #Continent
                alluvium=variable,fill=variable,label=variable,y=value)) + 
  scale_fill_manual(values = mycol,labels = custom_labels)+
  geom_stratum(aes(fill=variable),color='black',width = 0.5,size=0.1)+####width
  geom_flow(aes(fill=variable)) +
  guides(fill= guide_legend(ncol=1,reverse = F))+
  scale_y_continuous(limits = c(0,1.2),breaks = c(0,0.25,0.5,0.75,1)) +
  theme(panel.background = element_rect(fill='white', colour='white'), 
        panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x  = element_text(size=16, angle = 0,colour="black", face = "bold",hjust = 0.5),  
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 16, face = "bold"),
        strip.text=element_text(size=16))+
  labs(y = "Relative Abundance of COGs (%)", x="Group")+
  scale_y_continuous(expand = c(0,0))

ggsave('df3cog.pdf', width=10, height =6)



# se calculation following the previous mean calculation
df_se <- as.data.frame(t(apply(df, 1, function(x) {
  tapply(x, colnames(df), function(y) sd(y) / sqrt(length(y)))
})))

df_se <- df_se[, outord, drop = FALSE]

# Combine with means 
df_com <- cbind(cog = rownames(df1), mean = df1, se = df_se)
write.csv(df_com,'cog_mean_se.csv')


# COG circular barplot
# prepare input df with mean and SE

df<-read.csv('Features/cog_circular_dfa.csv')

ggplot(df, aes(x=as.factor(id), y=mean, fill=Var)) +       
  geom_bar(stat="identity", alpha=0.5, color="black", linewidth=0.1) +
  # change minus value to adjuste circle
  ylim(-5,12) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar()+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black", linewidth = 0.4)+
  # geom_text(aes(y = 0, label = COG), size = 3, angle = 0, color = "black")+
  scale_fill_manual(values = c("Clinical" = "#CAB2D6", 
                               "Environment" = "#7FC97F", 
                               "Food" = "#FBB4AE",
                               "Asia" = "#E5D8BD", 
                               "Europe" = "#DFC27D", 
                               "North_America" = "#BF812D",
                               "Food_Dry" = "#FB6A4A", 
                               "Food_Other" = "#E41A1C", 
                               "Food_Powdered" = "#FB9A99"))

ggsave('cog_circular_accessory.pdf', width=7, height=7)


### enriched VF, AMR genes, same method for VF and AMR
# dfa colnames as gene name
# prepare card_dfa file with dfa and card annotation file
card<-read.csv('card_dfa.csv', row.names=1, h=T)
  
meta <- read.csv("metadata/metadata.csv", h=T)
card$Bin<-rownames(card)
df1<-merge(card, meta, by ="Bin")

# revise df1 about sample/continent before next step
row.names(df1)<-df1[,1]
df1<-df1[,-1]
# ...

# make df1 numeric except for Continent and Sample_type
num <- sapply(df1, is.numeric)
num["Continent"] <- TRUE
# num["Sample_type"] <- TRUE
# num["Sample_type_specific"] <- TRUE

category_sums <- aggregate(. ~ Continent, df1[, num], sum)
#category_sums <- aggregate(. ~ Sample_type, df1[, num], sum)
#category_sums <- aggregate(. ~ Sample_type_specific, df1[, num], sum)

row.names(category_sums)<-category_sums[,1]
category_sums<-category_sums[,-1]
# remove South_America
category_sums<-category_sums[-4,]

# remove all 0 columns
all_zero_columns <- colSums(category_sums != 0) == 0
category_sums <- category_sums[, !all_zero_columns]


# prepare contingency tables for fisher test
contingency_tables <- list()

for (card in colnames(category_sums)) {
  contingency_table <- matrix(NA, nrow = 2, ncol = 3)
  rownames(contingency_table) <- c(card, paste0("Non-", card))
  colnames(contingency_table) <- rownames(category_sums)
  for (source in rownames(category_sums)) {
    contingency_table[1, source] <- category_sums[source, card]  
    contingency_table[2, source] <- sum(category_sums[source, ]) - category_sums[source, card]  
  }
  contingency_tables[[card]] <- contingency_table
}

## Fisher’s test
results <- data.frame(
  card = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each contingency table and perform Fisher's test with simulation
for (card in names(contingency_tables)) {
  fisher_result <- fisher.test(contingency_tables[[card]], simulate.p.value = TRUE, B = 1e5)
  results <- rbind(results, data.frame(
    card = card,
    P_value = fisher_result$p.value
  ))
}

results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")
write.csv(results,'card_results.csv')


### Random Forest
# dfa was generated by previous steps
dfa<-read.csv("dfaccessory.csv", row.names=1, h=T)
df1<-as.data.frame(t(dfa))

# Replace '~~~' (in the colname) with '_'
names(df1) <- gsub("~~~", "_", names(df1))  

sum(is.na(df1))

# check near zero variance and filter out, no metadata at this step
## default: freqCut = 95/5, uniqueCut = 10
nzv <- nearZeroVar(df1, saveMetrics = TRUE) 
df2 <- df1[, !nzv$nzv]

meta <- read.csv("metadata/metadata.csv")
df2$Bin <- rownames(df2)
df3<-merge(meta,df2, by="Bin")

#delete irrelated columns and assign row.names
row.names(df3)<-df3[,1]
df3 <- df3[, -1]
df4 <- df3[, -ncol(df3)]
# ...
# organize df accordingly and prepare for subsets 


# load organized subsets
# input files were generated using dfa and metadata
df<-read.csv('###.csv', row.names=1, h=T)
set.seed(123)
colnames(df)[1]<-c('group')
index <- createDataPartition(y = df$group, p = 0.7, list = FALSE)
train_set <- df[index, ]
test_set <- df[-index, ]

# data balancing
set.seed(123)
# Convert character columns to factors
train_set[sapply(train_set, is.character)] <- lapply(train_set[sapply(train_set, is.character)], as.factor) 
# Convert numeric columns to factors
train_set[sapply(train_set, is.numeric)] <- lapply(train_set[sapply(train_set, is.numeric)], as.factor)

table(train_set$group)

# 1. df Continent AEN
train_set.euro<-subset.data.frame(train_set, group=="Europe")
train_set2<-subset.data.frame(train_set, group!="Europe")
train_set2$group <- droplevels(train_set2$group)
train_smote <- SMOTE(group ~ ., data = train_set2, perc.over = 40, perc.under = 360)
train_smote2<-rbind(train_smote, train_set.euro)
# check for sample size
table(train_smote2$group)
table(train_set$group)
table(test_set$group)


# 2. df Sample_type CEF
train_smote <- SMOTE(group ~ ., data = train_set, perc.over = 100, perc.under = 400)

# 3. no smote for food subset

write.csv(train_set,'train_set.csv')
write.csv(test_set,'test_set.csv')
write.csv(train_smote,'train_smote.csv')

# random forest
# Random forest for classification expect the response variable to be a factor
# Categorical predictors should be factors (0 and 1 mean two categories)
# Numeric predictors should remain numeric

# Convert character columns to factors
train_set[sapply(train_set, is.character)] <- lapply(train_set[sapply(train_set, is.character)], as.factor) 
# Convert numeric columns to factors
train_set[sapply(train_set, is.numeric)] <- lapply(train_set[sapply(train_set, is.numeric)], as.factor)

# Convert for train_smote
# train_smote[sapply(train_smote, is.character)] <- lapply(train_smote[sapply(train_smote, is.character)], as.factor) 
# train_smote[sapply(train_smote, is.numeric)] <- lapply(train_smote[sapply(train_smote, is.numeric)], as.factor)

ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, search = "random", selectionFunction = "oneSE", classProbs = TRUE, allowParallel = TRUE)

# By default, caret uses 500 trees for the random forest model
rf_model <- train(Sample_type ~ ., data = train_set, method = "rf", trControl = ctrl, tuneLength = 10)

save(rf_model, file = "rf_model.RData")

## evaluate the model
# load rf model
rf_model<-rf_model_dfcAEN
print(rf_model)
# test set was generated by createDataPartition
test_set<-read.csv('test_set.csv')
# Convert character columns to factors
test_set[sapply(test_set, is.character)] <- lapply(test_set[sapply(test_set, is.character)], as.factor) 
# Convert numeric columns to factors
test_set[sapply(test_set, is.numeric)] <- lapply(test_set[sapply(test_set, is.numeric)], as.factor)

pred <- predict(rf_model, newdata = test_set)
conf_mat <- confusionMatrix(pred, test_set$Continent)
print(conf_mat)

# set scale = TRUE (scale 100); set scale = FALSE
df1 <- varImp(rf_model, scale = TRUE)	
df2 <- as.data.frame(df1$importance)
df2$Variable <- rownames(df2)
df3 <- df2[order(-df2$Overall),]
df4 <- head(df3, n = 50)
write.csv(df4,"top_50.csv")


# AUC curve
# Predict probabilities for each class instead of class predictions.
# Probabilities: Provide detailed insights into the model’s confidence and are essential for ROC curves, AUC, and threshold tuning.
# Class Predictions: Provide the final decision of the model for each instance and are used for calculating accuracy, confusion matrices, and deploying the model for practical applications.

predict <- predict(rf_model, newdata = test_set, type = "prob")
true_labels <- test_set$Continent

# Calculate ROC curves and AUC for each class One-vs-All Approach
roc_curves <- lapply(colnames(predict), function(class_name) {
  predict_class <- predict[, class_name]
  roc_curve <- roc(ifelse(true_labels == class_name, 1, 0), predict_class)
  auc_value <- auc(roc_curve)
  return(list(class_name = class_name, roc_curve = roc_curve, auc_value = auc_value))
})

for (i in seq_along(roc_curves)) {
  print(paste("AUC for", roc_curves[[i]]$class_name, ":", roc_curves[[i]]$auc_value))
}

# construct dataframe for figs
roc_data <- lapply(colnames(predict), function(class_name) {
  predict_class <- predict[, class_name]
  roc_curve <- roc(ifelse(true_labels == class_name, 1, 0), predict_class)
  data.frame(
    sensitivity = roc_curve$sensitivities,
    specificity = 1 - roc_curve$specificities,
    class = class_name
  )
})

# Combine ROC curve data
roc_data <- do.call(rbind, roc_data)

mycol <- c("Clinical" = "#CAB2D6", "Environment" = "#7FC97F", "Food" = "#FBB4AE")
mycol<-c("Asia" = "#E5D8BD", "Europe" = "#DFC27D", "North_America" = "#BF812D")
mycol <- c("Food_Dry" = "#FB6A4A", "Food_Other" = "#E41A1C", "Food_Powdered" = "#FB9A99")

roc_data$class <- factor(roc_data$class, levels = names(mycol))

custom_labels <- c("Clinical, auROC=###","Environment, auROC=###","Food, auROC=###")
custom_labels <- c("Asia, auROC=###","Europe, auROC=###","North America, auROC=###")
custom_labels <- c("Food_Dry, auROC=###", "Food_Other,auROC=###", "Food_Powdered,auROC=###")

ggplot(roc_data, aes(x = specificity, y = sensitivity, color = class)) +
  geom_line(linewidth = 1) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  theme(
    panel.background = element_blank(), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_line(color = "grey95"), 
    panel.border = element_rect(color = '#666666', fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 14), 
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(title = NULL)) +
  scale_color_manual(values = mycol,labels = custom_labels)

ggsave('df_ROC.pdf', width=7.5, height=4)

## top50 annotation
# remove "1" at the end of gene name
df4$Variable <- paste0(df4$Variable, "*")
df4$Variable <- sub("1\\*$", "", df4$Variable)

## prokka annotation
pan<-read.csv("Pan_analysis/pan_gene_annotation.csv")
colnames(df4)[2]<-c('Gene')
colnames(pan)[1]<-c('Gene')
df_pan<-merge(df4, pan, by="Gene")


## cog annotation
cog<-read.csv("Features/cog_annotation.csv")
colnames(cog)[1]<-c('Gene')
df_cog<-merge(df4, cog, by="Gene")


##VFDB annotation
hb<-read.csv("Features/vfdb_B_headers.csv")
sb<-read.csv("Features/vfdb_B_subject.csv")
colnames(sb)[1]<-c('Gene')
v1<-merge(df4,sb,by="Gene")
df_vf<-merge(v1,hb, by="subject")

##CARD
card<-read.csv("Features/card.csv")
colnames(card)[1]<-c('Gene')
df_card<-merge(df4, card, by="Gene")


