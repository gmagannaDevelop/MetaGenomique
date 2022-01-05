#!/usr/bin/env Rscript

rm(list=objects())
graphics.off()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quietly = TRUE)
BiocManager::install("qvalue")
list.of.packages <-c('devtools','RJSONIO','ecodist','gplots','scatterplot3d','usethis', 'httr','rcmdcheck','roxygen2','rversions')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(devtools)
library(tidyverse)
install_github(repo='MG-RAST/matR',quiet=T)
rm(list.of.packages,new.packages)

# Différents types de croûte : naturel, lavé, et fleuri.
# Un fromage est plus acide que le lait. On observe le PH.
# On a aussi du sel ajouté dans le fromage, qui va influencer la composition bactérienne.
# Pasteurisé => influence initialement les bactéries.
# On a des fromages de moutons, chèvres et vaches.

library(matR)

## Loading required package: MGRASTer
## MGRASTer (0.9 02e288)
## Loading required package: BIOM.utils
## BIOM.utils (0.9 dbcb27)
## matR: metagenomics analysis tools for R (0.9.1)
## List des accessions associées à l'étude mgp3362 (2 échantillons absents des métadonnées)

metadata.df <- read.table("TP1/MetadataCheese.csv", header = T)

list_mgp3362 <- metadata("mgp3362")$mgp3362
list_mgp3362 <- list_mgp3362[c(-11,-22)]

## Récupération des données taxonomiques (request='organism') depuis les hits
## de la base de données RDP (source='RDP') au niveau de l'ordre (group_level='order')
## avec une evalue de 1e-15 (evalue=15)
biom_phylum <- biomRequest(list_mgp3362, request="organism", source="RDP",
                           group_level="order", evalue=15, wait=TRUE)

exp.number <- ncol(biom_phylum)

clean_data_frame <- function(my.matrix){
  #' Convert a matrix to a dataframe, 
  #' eliminating null columns
  .df <- as.data.frame(my.matrix)
  .null.col.sums <- colSums(.df) == 0
  .to.drop <- .null.col.sums[.null.col.sums]
  .cols.to.drop <- names(.to.drop)
  return(.df[, !colnames(.df) %in% .cols.to.drop])
}

phylum.df <- clean_data_frame(as.matrix(biom_phylum))
phylum.probs <- apply(phylum.df, 2, function(x){ x / sum(x)})
OTU <- row.names(phylum.df)
names(OTU) <- NULL

# Calcul de la diversité alpha

compute_OTU_size <- function(phylum.exp){
  phylum.exp.probs <-phylum.exp/sum(phylum.exp)
  OTU.size <- c()
  sample.size <- seq(0,sum(phylum.exp), by=20)
  for (k in sample.size){
    OTU.size.k <- replicate(100, sample(OTU, k, prob = phylum.exp.probs, replace = TRUE) %>% unique() %>% length())
    OTU.size <- c(OTU.size,mean(OTU.size.k))
  }
  compute_OTU_size <- data.frame(sample.size, OTU.size)
}

rarefaction_curves <- list()

color.rarefaction_curves <- c(rep("blue",4),
                              "red",rep("blue",3),
                              rep("red",3),"blue",
                              rep("blue",3),"red",
                              "red",rep("blue",3),
                              "red","blue")

require(gridExtra)

for (k in 1:exp.number){
  
  OTU.size <- compute_OTU_size(phylum.df[,k])
  names(OTU.size) <- c("sample_size", "OTU_size")
  
  rarefaction_curves[[k]] <- ggplot(data = OTU.size, aes(x = sample_size, y = OTU_size)) +
    geom_line(color = color.rarefaction_curves[k]) +
    xlab(paste("Nombre de séquences (éch. ",k, ")", sep="")) + ylab("Nombre OTU") +
    theme_classic(base_size = 11)
  if (k%%4==0){
    grid.arrange(rarefaction_curves[[k-3]],
                 rarefaction_curves[[k-2]],
                 rarefaction_curves[[k-1]],
                 rarefaction_curves[[k]],
                 ncol=2, nrow=2,
                 top = "Diversité alpha")
  }
  if (k == exp.number){
    grid.arrange(rarefaction_curves[[k-1]],
                 rarefaction_curves[[k]],
                 ncol=2, nrow=2,
                 top = "Diversité alpha")
  }
}

rm(rarefaction_curves, color.rarefaction_curves)

# Calcul de l'entropie de Shannon

entropy <- function(probs.df){ 
  #' Shannon's entropy
  probs.clean <- probs.df[probs.df > 0] 
  - sum(probs.clean * log(probs.clean)) 
}

entropy <- apply(phylum.probs, 2, entropy)
print("The Shannon's entropies are:")
print(entropy)
rm(entropy)

## 2.2 Analyse des facteurs structurants la communauté

library(FactoMineR)
library(factoextra)

# phylum.df.normalize <- t(scale(as.matrix(phylum.df), scale = TRUE))
# res.pca <- PCA(phylum.df.normalize, scale.unit = FALSE, ncp = 6, graph = FALSE)

res.pca <- PCA(t(phylum.probs), scale.unit = FALSE, ncp = 4, graph = FALSE)

rind <- metadata.df$RindType
rind.df <- as.data.frame(rind)
rind.df <- as.factor(rind.df$rind)

fviz_pca_ind(res.pca,
             geom="point",
             axes = c(1, 2),
             pointsize=2,
             habillage = rind.df,
             invisible="quali",
) + labs(title="PCA (rind type)")

milk <- metadata.df$Milk
milk.df <- as.data.frame(milk)
milk.df <- as.factor(milk.df$milk)

fviz_pca_ind(res.pca,
             geom="point",
             axes = c(1, 2),
             pointsize=2,
             habillage = milk.df,
             invisible="quali",
) + labs(title="PCA (milk)")

pasteurized <- metadata.df$Pasteurized
pasteurized.df <- as.data.frame(pasteurized)
pasteurized.df <- as.factor(pasteurized.df$pasteurized)

fviz_pca_ind(res.pca,
             geom="point",
             axes = c(1, 2),
             pointsize=2,
             habillage = rind.df,
             invisible="quali",
) + labs(title="PCA (rind type)")

# LDA sans normalisation

library(MASS)

### Pasteurisation ###

res.lda.pasteurized <- lda(x=t(phylum.df), grouping = pasteurized)

decision.lda.pasteurized <- as.data.frame(colSums(res.lda.pasteurized$scaling * phylum.df))
names(decision.lda.pasteurized) <- "mean"

ggplot(data = decision.lda.pasteurized, aes(x = mean, fill=pasteurized)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.pasteurized$scaling[,1],names.arg = rownames(res.lda.pasteurized$scaling),las=2,
        col='black',cex.names = 0.76)

### Milk ###

res.lda.milk <- lda(x=t(phylum.df), grouping = milk)

decision.lda.milk <- as.data.frame(colSums(res.lda.milk$scaling * phylum.df))
names(decision.lda.milk) <- "mean"

ggplot(data = decision.lda.milk, aes(x = mean, fill=milk)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.milk$scaling[,1],names.arg = rownames(res.lda.milk$scaling),las=2,
        col='black',cex.names = 0.76)

### Milk ###

res.lda.rind <- lda(x=t(phylum.df), grouping = rind)

decision.lda.rind <- as.data.frame(colSums(res.lda.rind$scaling * phylum.df))
names(decision.lda.rind) <- "mean"

ggplot(data = decision.lda.rind, aes(x = mean, fill=rind)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Rind type")

par(mar=c(15,5,5,5))
barplot(res.lda.rind$scaling[,1],names.arg = rownames(res.lda.rind$scaling),las=2,
        col='black',cex.names = 0.76)

# LDA avec utilisation des quatre premières composantes principales

### Pasteurisation ###

res.lda.pasteurized <- lda(x=res.pca$ind$coord, grouping = pasteurized)

decision.lda.pasteurized <- as.data.frame(colSums(res.lda.pasteurized$scaling * phylum.df))
names(decision.lda.pasteurized) <- "mean"

ggplot(data = decision.lda.pasteurized, aes(x = mean, fill=pasteurized)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.pasteurized$scaling[,1],names.arg = rownames(res.lda.pasteurized$scaling),las=2,
        col='black',cex.names = 0.76)

### Milk ###

res.lda.milk <- lda(x=res.pca$ind$coord, grouping = milk)

decision.lda.milk <- as.data.frame(colSums(res.lda.milk$scaling * phylum.df))
names(decision.lda.milk) <- "mean"

ggplot(data = decision.lda.milk, aes(x = mean, fill=milk)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.milk$scaling[,1],names.arg = rownames(res.lda.milk$scaling),las=2,
        col='black',cex.names = 0.76)

### Rind ###

res.lda.rind <- lda(x=res.pca$ind$coord, grouping = rind)

decision.lda.rind <- as.data.frame(colSums(res.lda.rind$scaling * phylum.df))
names(decision.lda.rind) <- "mean"

ggplot(data = decision.lda.rind, aes(x = mean, fill=rind)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Rind type")

par(mar=c(15,5,5,5))
barplot(res.lda.rind$scaling[,1],names.arg = rownames(res.lda.rind$scaling),las=2,
        col='black',cex.names = 0.76)


# LDA avec normalisation

### Pasteurisation ###

res.lda.rind <- lda(x=t(phylum.df), grouping = rind)

phylum.normalize.df <- t(scale(as.matrix(phylum.df), scale = TRUE))

res.lda.pasteurized <- lda(x=t(phylum.normalize.df), grouping = pasteurized)

decision.lda.pasteurized <- as.data.frame(colSums(res.lda.pasteurized$scaling * phylum.df))
names(decision.lda.pasteurized) <- "mean"

ggplot(data = decision.lda.pasteurized, aes(x = mean, fill=pasteurized)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_classic(base_size = 11) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.pasteurized$scaling[,1],names.arg = rownames(res.lda.pasteurized$scaling),las=2,
        col='black',cex.names = 0.76)

### Milk ###

res.lda.milk <- lda(x=res.pca$ind$coord, grouping = milk)

decision.lda.milk <- as.data.frame(colSums(res.lda.milk$scaling * phylum.df))
names(decision.lda.milk) <- "mean"

ggplot(data = decision.lda.milk, aes(x = mean, fill=milk)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Pasteurized")

par(mar=c(15,5,5,5))
barplot(res.lda.milk$scaling[,1],names.arg = rownames(res.lda.milk$scaling),las=2,
        col='black',cex.names = 0.76)

### Rind ###

res.lda.rind <- lda(x=res.pca$ind$coord, grouping = rind)

decision.lda.rind <- as.data.frame(colSums(res.lda.rind$scaling * phylum.df))
names(decision.lda.rind) <- "mean"

ggplot(data = decision.lda.rind, aes(x = mean, fill=rind)) + 
  geom_histogram(aes(y=..density..), color="#e9ecef", alpha=0.6, position = 'identity', bins = 18) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "red")) +
  theme_classic(base_size = 15) +
  labs(fill="Rind type")

par(mar=c(15,5,5,5))
barplot(res.lda.rind$scaling[,1],names.arg = rownames(res.lda.rind$scaling),las=2,
        col='black',cex.names = 0.76)

