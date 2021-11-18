## ----setup, include=FALSE------------------------------------------------------------------------
knitr::opts_chunk$set(echo = F, warning = F, message = F)


## ----lib.import----------------------------------------------------------------------------------
library(matR)
library(tidyverse)
library(reshape2)


## ----data.import---------------------------------------------------------------------------------
metadata.df <- read.table("MetadataCheese.csv", header = T)
metadata.tb <- tibble(metadata.df)

#List des accessions associées à l'étude mgp3362 
# (2 échantillons absents des métadonnées)
list_mgp3362<-metadata("mgp3362")$mgp3362
list_mgp3362<-list_mgp3362[c(-11,-22)]
#Récupération des données taxonomiques (request='organism') depuis les hits
#de la base de données RDP (source='RDP') au niveau de l'ordre (group_level='order')
#avec une evalue de 1e-15 (evalue=15)
biom_phylum<-biomRequest(list_mgp3362,request="organism",source="RDP",
group_level="order",evalue=15,wait=TRUE)

#Transformation en matrice
phylum.matrix <- as.matrix(biom_phylum)
clean_data_frame <- function(my.matrix){
  #' Convert a matrix to a dataframe, 
  #' eliminating null columns
  .df <- as.data.frame(my.matrix)
  .null.col.sums <- colSums(.df) == 0
  .to.drop <- .null.col.sums[.null.col.sums]
  .cols.to.drop <- names(.to.drop)
  return(.df[, !colnames(.df) %in% .cols.to.drop])
}

phylum.df <- clean_data_frame(phylum.matrix)
phylum.tb <- phylum.df %>% rownames_to_column("Strain") %>% as_tibble()


## ----rarefaction.curve---------------------------------------------------------------------------
otus <- phylum.tb %>% select(Strain) %>% unlist
names(otus) <- NULL
phylum.prob.df <- apply(phylum.df, 2, function(x){ x / sum(x)})


rarefactions.ls <- list()
for (i in 1:ncol(phylum.df)){
  .n.taxa <- c()
  .n.seq <- c()
  for (j in seq(1, 10000, 50)){
    .tax <- sample(
      otus, j, prob = phylum.prob.df[,i], replace = T
    ) %>% unique() %>% length()
    .n.taxa <- c(.n.taxa, .tax)
    .n.seq <- c(.n.seq, j)
  }
  rarefactions.ls[[i]] <- tibble(nSeq=.n.seq, nOTU=.n.taxa)
}

rarefactions.ls[[2]] %>% 
  ggplot(aes(x=nSeq, y=nOTU)) + geom_line() + geom_smooth()


## ----shannon-------------------------------------------------------------------------------------
H <- function(x){ 
  #' Shannon's entropy
  p <- x[x > 0] 
  - sum(p * log(p)) 
}

apply(phylum.prob.df, 2, H)


## ----pca.lda-------------------------------------------------------------------------------------
#Chargement du package MASS
library(MASS)
## NB la matrice de données doit être transformée
LDA<-lda(x=t(phylum.matrix),grouping=metadata.df$Pasteurized)

## Calcul des valeurs pour chaque groupe
LDA1_RawMilk <- colSums(apply(phylum.matrix[,metadata.df$Pasteurized=='N'],2,
function(x){LDA$scaling*x}))

LDA1_Pasteurized<-colSums(apply(phylum.matrix[,metadata.df$Pasteurized=='Y'],2,
function(x){LDA$scaling*x}))

#Représentation sous forme d'histogramme
plot.new()
hist(LDA1_Pasteurized, xlim=c(min(LDA1_RawMilk)-4, max(LDA1_Pasteurized)+1),
     col='green', xlab='LD1 means', ylab='Frequency', main='Predicted Values'
)

hist(LDA1_RawMilk,col='red',add=T)
legend('topleft',pch=19,col=c('red','green'),legend = c('Raw','Pasteurized'))



## ----other.plot----------------------------------------------------------------------------------
par(mar=c(15,5,5,5))
barplot(LDA$scaling[,1],names.arg = rownames(LDA$scaling),las=2,
col='black',cex.names = 0.76)


## ----acp-----------------------------------------------------------------------------------------
#lda(x=t(phylum_matrix),grouping=metadata$Pasteurized)
pca_res <- acp <- prcomp(t(phylum.prob.df), scale = F)
biplot(acp)
acp$x[,1:2] %>% as.data.frame() %>% ggplot(aes(PC1, PC2)) + geom_point()

var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
pca_res$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(size=4) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

