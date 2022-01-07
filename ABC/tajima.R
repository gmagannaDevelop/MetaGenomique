#!/usr/bin/env Rscript

#' ---
#' title: "ABC"
#' author:
#'   - Gustavo MAGAÑA LOPEZ
#'   - Théo RONCALLI
#' output:
#'   png_document
#' ---
#'

#### Loading dependencies ####

rm(list=objects())
graphics.off()

.dependencies <- 
  c("here", "abc", "coala", 
    "magrittr", "ggplot2", "dplyr",
    "tibble", "foreach", "doRNG",
    "RColorBrewer", "ggplot2"
  )
new.packages <- 
  .dependencies[!(.dependencies %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)
.. <- sapply(.dependencies, function(x){ library(x, character.only = T)})

source(here("ABC/Fonctions_TP_demo_real.R"))

#### Loading Data ####

hap.fname <- here("ABC/chr22.CG_54genomes_shapeit_phased.txt.haps")
info.fname <- here("ABC/CG_54genomes_indiv.txt")
realdata <- read.hap.geno(hap.input=hap.fname, info.input = info.fname)


names(realdata)
realdata[["info"]][["POP"]] %>% unique()  # Toutes les populations recensées dans notre jeu de données

#### Color initialization ####

code.population <- unique(realdata$info[,"POP"])
code.region <- unique(realdata$info[,"REGION"])

color.population <- brewer.pal(length(code.population), "Paired")
names(color.population) <- code.population

color.region <- brewer.pal(length(code.region), "Set1")
names(color.region) <- code.region

#### Statistics computation ####

stat.1000g <- NULL
for (pop in unique(realdata$info$POP)) {
  # subset the dataset
  pop.subset <- select.label(realdata, label = pop)
  # compute estimates
  
  size.norm.factor <- pop.subset$region_bp_lim[2] - pop.subset$region_bp_lim[1]
  theta_watterson <- calculate.thetaW(pop.subset$segsites) / size.norm.factor
  theta_tajima <- calculate.thetaPi(pop.subset$segsites) / size.norm.factor
  tajimas_d <- calculate.Dt(pop.subset$segsites)
  tajimas_d_var <- calculate.Dt.var(pop.subset$segsites)
  # create the dataframe
  .pop.stats <- 
    data.frame(
      Theta.Tajima=theta_tajima, D.Tajima=tajimas_d, 
      Var.D.Tajima=tajimas_d_var,
      Theta.Watterson=theta_watterson
    )
  stat.1000g <- rbind(stat.1000g, .pop.stats)
}

rownames(stat.1000g) <- unique(realdata$info$POP)

sumstats <- stat.1000g %>% 
  rownames_to_column("Population") %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(color = color.population)

plot.theta_tajima <- ggplot(sumstats, aes(x = Population, y = Theta.Tajima, fill = color)) +
  geom_bar(stat = "identity") +
  ylab("Theta Tajima") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

paste(here("ABC/Figures"), "Theta_Tajima.png", sep="/") %>%
  ggsave(plot = plot.theta_tajima, width = 6, height = 5)
rm(plot.theta_tajima)

plot.theta_watterson <- ggplot(sumstats, aes(x = Population, y = Theta.Watterson, fill = color)) +
  geom_bar(stat = "identity") +
  ylab("Theta Watterson") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

paste(here("ABC/Figures"), "Theta_Watterson.png", sep="/") %>%
  ggsave(plot = plot.theta_watterson, width = 6, height = 5)
rm(plot.theta_watterson)

plot.D_tajima <- ggplot(sumstats, aes(x = Population, y = D.Tajima, fill = color)) +
  geom_bar(stat = "identity") +
  ylab("D Tajima") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

paste(here("ABC/Figures"), "D_Tajima.png", sep="/") %>%
  ggsave(plot = plot.D_tajima, width = 6, height = 5)
rm(plot.D_tajima)
