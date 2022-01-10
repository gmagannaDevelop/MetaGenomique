#' ---
#' title: "ABC"
#' author:
#'   - Gustavo MAGAÑA LOPEZ
#'   - Théo RONCALLI
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---
#'

#### Loading dependencies ####
{
  # declare dependencies
  .dependencies <- 
    c("here", "abc", "coala", "reshape2",
      "magrittr", "ggplot2", "dplyr", "grid",
      "tibble", "foreach", "doRNG", "stringr",
      "RColorBrewer", "gridExtra"
      )
  # Check if some are not installed
  new.packages <- 
    .dependencies[!(.dependencies %in% installed.packages()[,"Package"])]
  # Install missing packages
  if(length(new.packages)) install.packages(new.packages)
  # Silently load all the dependencies
  .. <- sapply(.dependencies, function(x){ library(x, character.only = T)})
  
  # load utility functions
  source(here("ABC/Fonctions_TP_demo_real.R")) # this script has slighly changed since TP1
  source(here("ABC/customfuncs.R")) # To compute tables, error metrics, etc
  # TODO : merge cv4postpr and cv4abc into a single parallel_abc file 
  source(here("ABC/cv4postpr.R")) # To access parallel implementations of the functions
}

# load data
{
  hap.fname <- here("ABC/chr22.CG_54genomes_shapeit_phased.txt.haps")
  info.fname <- here("ABC/CG_54genomes_indiv.txt")
  realdata <- read.hap.geno(hap.input=hap.fname, info.input = info.fname)
}
set.seed(1234)


#' 2.2 Calcul de statistiques résumées par population.
{
  stat.1000g.df <- NULL
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
    stat.1000g.df <- rbind(stat.1000g.df, .pop.stats)
 }
  rownames(stat.1000g.df) <- unique(realdata$info$POP)
  # We realised that simulated data only has three statistics :
  # (Pi, D Tajima et Var(D Tajima)), so we take the subset
  stat.1000g <- stat.1000g.df[, 1:3]
  colnames(stat.1000g) <- c("pi", "tajimas_d", "tajimas_d_var")
  # the tibble facilitates using tidyverse's tools
  stat.1001g.tb <- stat.1000g %>% rownames_to_column("Population") %>% tibble()
}

#' Simulated data
loaded <- load(file = here("ABC/simu.expeA.RData"))
loaded
# [1] "sumstats.expeA" "models.expeA"   "params.bott"
# Three objects were imported from the file

#' Data analysis
dim(sumstats.expeA)
head(sumstats.expeA) # each line is a set of parameters... for a different model ?
table(models.expeA) # equally divided bott, const, exp -> 50_000 each
head(models.expeA)
dim(params.bott)  
head(params.bott)

# Visualise model parameters
{
aggregated.expeA <- cbind(sumstats.expeA, models.expeA)
colnames(aggregated.expeA) <- c("Theta.Tajima", "D.Tajima", "Var.D.Tajima", "Type")

boxplot.theta_tajima <- ggplot(aggregated.expeA, aes(x=Type, y=Theta.Tajima)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE,
               fill=brewer.pal(3, "Dark2"), color="black") +
  ylab("Theta Tajima") +
  xlab("Type") +
  theme_classic()

paste(here("ABC/Figures"), "boxplot_Theta_Tajima.png", sep="/") %>%
  ggsave(plot = boxplot.theta_tajima, width = 6, height = 5)

boxplot.d_tajima <- ggplot(aggregated.expeA, aes(x=Type, y=D.Tajima)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE,
               fill=brewer.pal(3, "Dark2"), color="black") +
  ylab("D Tajima") +
  xlab("Type") +
  theme_classic()

paste(here("ABC/Figures"), "boxplot_D_Tajima.png", sep="/") %>%
  ggsave(plot = boxplot.d_tajima, width = 6, height = 5)

boxplot.var_d_tajima <- ggplot(aggregated.expeA, aes(x=Type, y=Var.D.Tajima)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE,
               fill=brewer.pal(3, "Dark2"), color="black") +
  ylab("Variance (D Tajima)") +
  xlab("Type") +
  theme_classic()

paste(here("ABC/Figures"), "boxplot_Var_D_Tajima.png", sep="/") %>%
  ggsave(plot = boxplot.var_d_tajima, width = 6, height = 5)

aggregated_boxplot_tajima <- grid.arrange(boxplot.theta_tajima, boxplot.d_tajima, boxplot.var_d_tajima, ncol=3)

paste(here("ABC/Figures"), "Aggregated_boxplot_tajima.png", sep="/") %>%
  ggsave(plot = aggregated_boxplot_tajima, width = 6, height = 5)

rm(boxplot.theta_tajima, boxplot.d_tajima, boxplot.var_d_tajima, aggregated_boxplot_tajima)
}


# Explore main ABC function
#?postpr

# Explore how the functions should be used, as well as its outputs.
{
  random_idx <- 
    sample.int(length(models.expeA), size = 1)
  random_idx
  
  abc.mod.sel <- postpr(
    target = sumstats.expeA[random_idx,],   # pseudo-observed sumstats
    index = models.expeA[-random_idx],      # known models for reference db
    sumstat = sumstats.expeA[-random_idx,], # sumstats for reference db
    method = 'rejection',                   # type of ABC algorithm
    tol=0.05                                # percentage of kept simulations
  )
  prob <- summary(abc.mod.sel)$Prob
  prob
}

#' 3.3 Sélection de modèle démographique, évaluation des performances
#' This shoult be performed using the script we developped for this purpose
#' From the repo's root, you should be able to run it as follows :
#' `Rscript ABC/parallel_cv4postpr.R 150 0.005 0.1 20 14`
#' (if called with no arguments it will automatically print a help message)

#' If you really want to run it on your R session, the function should
#' already be loaded into the global namespace, you can call it as follows :
#' 
# parallel.cv.modsel.reject <- 
#   parallel_cv4postpr(
#     models.expeA, sumstats.expeA, nval=5, tols=seq(0.005, 0.95, length.out = 20)
#   )

#' Gridsearch
{
  # Load data from parallel_cv4postpr.R runs
  general <- my.read.table(here("ABC/data/cv4postpr_nval_150_tols_0.005_0.95_50.csv"))
  specific <- my.read.table(here("ABC/data/cv4postpr_nval_150_tols_0.005_0.1_50.csv"))
  
  # Visualise the results
  general.plot <- plot_cv_rej(general) + geom_smooth(method="lm")
  specific.plot <- plot_cv_rej(specific) + geom_smooth()
}

#' Calculate goodnes-of-fit
{
  #estimCEU <- stat.1000g[population, ]
  pop.names <- as.list(rownames(stat.1000g)) 
  names(pop.names) <- pop.names
  gofs <- 
    sapply(pop.names, gof_for_population, USE.NAMES = T, simplify = F)
}

likeliest.scenario <- data.frame(
   Population = names(gofs),
   P.bott = as.numeric(sapply(gofs, function(x){ as.numeric(x$posterior.probs['bott']) })),
   P.const = as.numeric(sapply(gofs, function(x){ as.numeric(x$posterior.probs['const']) })),
   P.exp = as.numeric(sapply(gofs, function(x){ as.numeric(x$posterior.probs['exp']) })),
   Scenario = as.character(sapply(gofs, 
     function(x){ 
       names(which.max(x$posterior.probs)) 
     }, 
     USE.NAMES = T, 
     simplify = F
   ))
)
likeliest.scenario
my.write.csv(likeliest.scenario, filename = here("ABC/data/likeliest.csv"))

#' # 4. Estimer les paramètres d’un modèle donné
sumstats.bott <- subset(sumstats.expeA, subset=models.expeA=="bott")
colnames(params.bott)
params.bott.tb <- tibble(params.bott)

suppressMessages({
  .ne.hist <- params.bott.tb %>% ggplot(aes(x=Ne)) + geom_histogram(colour="black", fill="white")
  .a.hist <- params.bott.tb %>% ggplot(aes(x=a)) + geom_histogram(colour="black", fill="white")
  .dur.hist <- params.bott.tb %>% ggplot(aes(x=duration)) + geom_histogram(colour="black", fill="white")
  .star.hist <- params.bott.tb %>% ggplot(aes(x=start)) + geom_histogram(colour="black", fill="white")
  .param.hists <- grid.arrange(
    .ne.hist, .a.hist, .dur.hist, .star.hist, 
    nrow=2, ncol=2,
    top=textGrob("Distribution empirique des paramètres")
  )
  .param.hists
})

{
sumstats.bott.tb <- tibble(sumstats.bott) %>% 
  dplyr::select(
    `Theta Tajima (pi)` = pi,
    `D Tajima` = tajimas_d,
    `Var(D Tajima)` = tajimas_d_var
  )
sumstats.bott.tb$Ne <- params.bott$Ne

  .pi.plot <- sumstats.bott.tb %>% 
    ggplot(aes(x=Ne, y=`Theta Tajima (pi)`)) + geom_point(alpha=0.03) + geom_smooth(alpha=0.5)
  .d.plot <- sumstats.bott.tb %>% 
    ggplot(aes(x=Ne, y=`D Tajima`)) + geom_point(alpha=0.03) + geom_smooth(alpha=0.5)
  .var.d.plot <- sumstats.bott.tb %>% 
    ggplot(aes(x=Ne, y=`Var(D Tajima)`)) + geom_point(alpha=0.03) + geom_smooth(alpha=0.5)
  .param.trends <- gridExtra::grid.arrange(
    .pi.plot, .d.plot, .var.d.plot, 
    nrow=3, ncol=1,
    top=textGrob("Statistiques ~ Ne")
  )
  .param.trends
}
# Tracez les statistiques résumées en fonction du paramètre de taille efficace Ne, 
# commentez. Pourquoi ces statistiques résumées vous paraissent utiles pour 
# réaliser l’étude ABC complète à savoir 
# sélection de modèle + inférence de paramètres.

#' # 4.1 Exploration, intuitions
#' # 4.2 Exemple d’inférence
# ?abc
{
  random_idx <- sample.int(length(models.expeA[models.expeA=="bott"]), size = 1)
  random_idx
  res <- abc(target = sumstats.bott[random_idx,],
             param = params.bott,  
             sumstat =sumstats.bott[-random_idx],
             tol=0.05,
             method = "loclinear")
  
  cat(dim(res$unadj.values), "\n")
  cat(head(res$unadj.values), "\n")
  cat(summary(res$pi), "\n")
  plot_abc_distributions(res, params.bott, "Ne")
}

#' This shoult be performed using the script we developped for this purpose
#' From the repo's root, you should be able to run it as follows :
#' `Rscript ABC/parallel_cv4abc.R 150 0.005 0.1 20 16 ridge` (takes about an hour)
#' (if called with no arguments it will automatically print a help message)

#' If you really want to run it on your R session, the function should
#' already be loaded into the global namespace, you can call it as follows :
#' 
# # sequential, for exploration
# cv <- cv4abc(
#   param = params.bott, sumstat = sumstats.bott,
#   nval=10, tol=c(0.05, 0.1), method="loclinear", abc.out = NULL
# )
# 
# # the real thing ;)
# cv2 <- parallel_cv4abc(param = params.bott, sumstat = sumstats.bott,
#                 nval=10, tols = c(0.05, 0.1),
#                 method = "loclinear", abc.out = NULL, nthreads = 12)

global.abc <- my.read.table("ABC/data/cv4abc_loclinear_nval_150_tols_0.005_0.95_50.csv")
glob.plot <- plot_abc_cv_rmse(global.abc) + labs(title = "Global exploration")

# Gridsearch
{
  # Load data from parallel_cv4abc.R runs
  cv.loclinear <- 
    my.read.table(here("ABC/data/cv4abc_loclinear_nval_150_tols_0.005_0.3_20.csv"))
  cv.ridge <- 
    my.read.table(here("ABC/data/cv4abc_ridge_nval_150_tols_0.005_0.3_20.csv"))
  
  # Visualise the results
  abc.loclinear.cv.plot <- plot_abc_cv_rmse(cv.loclinear) + labs(title="method = 'loclinear'")
  abc.ridge.cv.plot <- plot_abc_cv_rmse(cv.ridge) + labs(title="method = 'ridge'")
  narrow.abc <- grid.arrange(
    abc.loclinear.cv.plot, abc.ridge.cv.plot, 
    nrow=2, ncol=1,
    top=textGrob("RMSE across parameters")
  )
}
# how to account for different behaviours ? 
# compute a weighted mean !

{
# caculated normalised RMSE (scale(center=F, scale=T))
cv.loclinear.norm <- normalise_rmse(cv.loclinear, .melt = F)
cv.ridge.norm <- normalise_rmse(cv.ridge, .melt = F)
# compute weighted averages (using their standard deviation)
cv.loclinear.w <- rmse_weighted_average(cv.loclinear.norm)
cv.ridge.w <- rmse_weighted_average(cv.ridge.norm)
# compare methods
cv.methods <- data.frame(
  tolerance=cv.loclinear.w$Tolerance,
  loclinear=cv.loclinear.w$w.RMSE,
  ridge=cv.ridge.w$w.RMSE
) %>% 
  reshape2::melt(
    id.vars="tolerance", 
    variable.name = "Method",
    value.name = "Scaled Error"
  )

method.comparison <- cv.methods %>% 
  ggplot(aes(x=tolerance, y=`Scaled Error`)) + 
  geom_point(aes(colour=Method)) + 
  labs(title="Weighted average of RMSE")
}


#' # 4.4 Application aux données réelles
#' 
#' 

# Select bottleneck populations
{
pops.bott <- 
  likeliest.scenario[ likeliest.scenario$Scenario == "bott", "Population"]
names(pops.bott) <- pops.bott
}


test1 <- sapply(
    pops.bott, 
    abc_with_distributions,
    col.name="Ne",
    method="ridge",
    tol = 0.025,
    USE.NAMES = T, 
    simplify = F
)

{
the.params <- colnames(params.bott)
names(the.params) <- the.params
}

real.abc <- sapply(the.params, function(x){
  sapply(pops.bott, abc_with_distributions,
  col.name=x, method="ridge", tol = 0.025,
  USE.NAMES = T, simplify = F)
  }, USE.NAMES = T, simplify = F
)

grid.arrange(real.abc$Ne$CEU, real.abc$a$CEU, real.abc$duration$CEU, real.abc$start$CEU, nrow=2, ncol=2)

{
resCHB <- abc(target = stat.1000g["CHB", ],
              param = params.bott,  
              sumstat =sumstats.bott,
              tol=0.025,
              method = "ridge")
  
chb.res <- list()
chb.res[[1]] <- plot_abc_distributions(res.abc = resCHB, params.df = params.bott, col.name = "Ne")#abc_with_distributions("CHB", "Ne", tol = 0.025, method = "ridge")
chb.res[[2]] <- plot_abc_distributions(res.abc = resCHB, params.df = params.bott, col.name = "a")#abc_with_distributions("CHB", "a", tol = 0.025, method = "ridge")
chb.res[[3]] <- plot_abc_distributions(res.abc = resCHB, params.df = params.bott, col.name = "duration")#abc_with_distributions("CHB", "duration", tol = 0.025, method = "ridge")
chb.res[[4]] <- plot_abc_distributions(res.abc = resCHB, params.df = params.bott, col.name = "start")#abc_with_distributions("CHB", "start", tol = 0.025, method = "ridge")
plotCHB <- do.call(grid.arrange, chb.res)
}

{
  .estim.ss <- resCHB$ss %>% apply(2, mean)
  .obs.ss <- stat.1000g["CHB", ]
  chb.df <- rbind(.estim.ss, .obs.ss)
  rownames(chb.df) <- c("Estimate", "Observed")
  chb.df %>% rownames_to_column("Type") %>% my.write.csv(filename = here("ABC/data/chb.csv"))
}

