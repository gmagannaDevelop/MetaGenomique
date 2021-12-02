#Simulation forward drift-mutation
Simul.forward<-function(Nsim,Ngen,PopSize,mu,Finit){
  #Nsim : number of simulations/observations, replicates
  #Ngen : number of generations
  #PopSize, Ne
  #mu, mutation rate
  #Finit, initial frequency
  Pop<-rep(Finit,Nsim)
  Pop<-rbind(Pop,rbinom(Nsim,PopSize,Pop)/PopSize)
  Mutation<-rep(sample(c(-1,1),1),Nsim)*rpois(Nsim,PopSize*mu)/PopSize
  Pop[2,]<-as.numeric(Pop[2,])+Mutation
  Pop[2,]<-ifelse(Pop[2,]>1,1,ifelse(Pop[2,]<0,0,Pop[2,]))
  gen<-2
  while(gen<Ngen+1){
    Pop<-rbind(Pop,rbinom(Nsim,PopSize,Pop[gen,])/PopSize)			
    Mutation<-rep(sample(c(-1,1),1),Nsim)*rpois(Nsim,PopSize*mu)/PopSize
    Pop[length(Pop[,1]),]<-Pop[length(Pop[,1]),]+Mutation
    Pop[length(Pop[,1]),]<-ifelse(Pop[length(Pop[,1]),]>1,1,ifelse(Pop[length(Pop[,1]),]<0,0,Pop[length(Pop[,1]),]))
    gen<-gen+1
  }
  plot(1:(Ngen+1),Pop[,1],type="l",xlab="Generations",ylab="Allele frequency",ylim=c(0,1))
  apply(Pop[,2:Nsim],2,function(Z){lines(1:(Ngen+1),Z,col=sample(rainbow(1000),1))})
}


##Coalescence simulations
#loading msms in coala
list.of.packages <- c("ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages)
}
#Download these packages from ecampus and install as follow
#install.packages('scrm_1.7.3-1.tar.gz',repos=NULL,type='source')
#install.packages('assertthat_0.2.1.tar.gz',repos=NULL,type='source')
#install.packages("rehh")
#install.packages("RcppArmadillo")
#install.packages('coala_0.6.0.tar.gz',repos=NULL,type='source')
library(coala) #for coalescent models
library(ape) # for coalescent tree reconstruction

activate_msms("msms3.2rc-b163.jar",priority = 1000)

#neutral model
Nlines=10 #Number of lineages ie haplotypes
K=100 #Number of chromosomes ie independant loci
Theta=10 #Theta = 4Neu
model_neut <- coal_model(Nlines,K) +
  feat_mutation(Theta) +
  sumstat_sfs() +
  sumstat_tajimas_d() +
  sumstat_nucleotide_div() +
  sumstat_seg_sites()+
  sumstat_trees()

#launch model
stats_neut<-simulate(model_neut)
#Watterson's theta
watterson_neut<-unlist(lapply(stats_neut$seg_sites,function(x){dim(x)[2]/sum(1/1:(sum(Nlines)-1))}))

#plot output
par(mfrow=c(3,2))
hist(stats_neut$pi,main="Tajima's Theta")
hist(watterson_neut,main="Watterson's Theta")
hist(stats_neut$tajimas_d,main="Tajima's D")
barplot(stats_neut$sfs,main="Derived site frequency spectrum")

#plot 1st coalescent tree
plot.phylo(read.tree(text=stats_neut$trees[[1]]))



#Exponential growth model
Nlines=10 #Number of lineages ie haplotypes
K=100 #Number of chromosomes ie independant loci
Theta=10 #Theta = 4Neu
alpha= 10 #Exponential growth (ie exponential decrease backward ) of the population following Nt=No*exp(-alpha*t) with Nt = effective population size back at time t and No current effective population size
model_growth <- coal_model(Nlines,K) +
  feat_mutation(Theta) +
  feat_growth(alpha)+
  sumstat_sfs() +
  sumstat_tajimas_d() +
  sumstat_nucleotide_div() +
  sumstat_seg_sites()+
  sumstat_trees()

#launch model
stats_growth<-simulate(model_growth)
#Watterson's theta
watterson_growth<-unlist(lapply(stats_growth$seg_sites,function(x){dim(x)[2]/sum(1/1:(sum(Nlines)-1))}))

#plot output
par(mfrow=c(3,2))
hist(stats_growth$pi,main="Tajima's Theta")
hist(watterson_growth,main="Watterson's Theta")
hist(stats_growth$tajimas_d,main="Tajima's D")
barplot(stats_growth$sfs,main="Derived site frequency spectrum")

#plot 1st coalescent tree
plot.phylo(read.tree(text=stats_growth$trees[[1]]))
plot.phylo(read.tree(text=stats_growth$trees[[2]]))


#Bottleneck model
Nlines=10 #Number of lineages ie haplotypes
K=100 #Number of chromosomes ie independant loci
Theta=10 #Theta = 4Neu
fraction=0.01 #Strength of the bottleneck as a proportion of the ancestral effective population size ie 0.5 means current size is 50% of ancestral popsize
model_bottleneck <- coal_model(Nlines,K) +
  feat_mutation(Theta) +
  feat_size_change(fraction,time=0)+
  sumstat_sfs() +
  sumstat_tajimas_d() +
  sumstat_nucleotide_div() +
  sumstat_seg_sites() +
  sumstat_trees()

#launch model
stats_bottleneck<-simulate(model_bottleneck)
#Watterson's theta
watterson_bottleneck<-unlist(lapply(stats_bottleneck$seg_sites,function(x){dim(x)[2]/sum(1/1:(sum(Nlines)-1))}))
#plot output
par(mfrow=c(2,2))
hist(stats_bottleneck$pi,main="Tajima's Theta")
hist(watterson_bottleneck,main="Watterson's Theta")
hist(stats_bottleneck$tajimas_d,main="Tajima's D")
barplot(stats_bottleneck$sfs,main="Derived site frequency spectrum")

#plot 1st coalescent tree
dev.off();plot.phylo(read.tree(text=stats_bottleneck$trees[[1]]))

#Two populations split model
Nlines=c(10,20) #Number of lineages ie haplotypes in each populations
K=100 #Number of chromosomes ie independant loci
Theta=100 #Theta = 4Neu
merge_time=5 #Time back in time when population marge as a fraction of 2Ne
source_pop=1 # The source population
splitin_pop=2 # The population that split from the source population
model_split <- coal_model(Nlines,K) +
  feat_mutation(Theta) +
  feat_pop_merge(merge_time,source_pop,splitin_pop)+
  sumstat_sfs() +
  sumstat_tajimas_d() +
  sumstat_nucleotide_div() +
  sumstat_seg_sites() +
  sumstat_trees()

#launch model
stats_split<-simulate(model_split)
#Watterson's theta -> problem overestimated!!!
watterson_split<-(unlist(lapply(stats_split$seg_sites,function(x){(dim(x)[2]-sum(colSums(as.matrix(x))==Nlines[1] | colSums(as.matrix(x))==Nlines[2]))/sum(1/1:(sum(Nlines)-1))}))/2)

#plot output
par(mfrow=c(2,2))
hist(stats_split$pi,main="Tajima's Theta")
hist(watterson_split,main="Watterson's Theta")
hist(stats_split$tajimas_d,main="Tajima's D")
barplot(stats_split$sfs,main="Derived site frequency spectrum")

#plot 1st coalescent tree
dev.off();plot.phylo(read.tree(text=stats_split$trees[[1]]))

#neutral model with recombination
Nlines=10 #Number of lineages ie haplotypes
K=100 #Number of chromosomes ie independant loci
Theta=10 #Theta = 4Neu
Rho=1000 #Rho=4ner
model_neutrec <- coal_model(Nlines,K) +
  feat_mutation(Theta) +
  feat_recombination(Rho) +
  sumstat_sfs() +
  sumstat_tajimas_d() +
  sumstat_nucleotide_div() +
  sumstat_seg_sites()+
  sumstat_trees()

#launch model
stats_neutrec<-simulate(model_neutrec)
#Watterson's theta
watterson_neutrec<-unlist(lapply(stats_neutrec$seg_sites,function(x){dim(x)[2]/sum(1/1:(sum(Nlines)-1))}))

#plot output
par(mfrow=c(2,2))
hist(stats_neutrec$pi,main="Tajima's Theta")
hist(watterson_neutrec,main="Watterson's Theta")
hist(stats_neutrec$tajimas_d,main="Tajima's D")
barplot(stats_neutrec$sfs,main="Derived site frequency spectrum")

#plot 1st coalescent tree
dev.off();plot.phylo(read.tree(text=stats_neutrec$trees[[1]][1]))

