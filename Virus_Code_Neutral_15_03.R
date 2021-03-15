#---
#  title: "Virus replication simulation"
#author: "Amicone Massimo"
#date: "10/2/2021"
#output: pdf_document
#
#---
  
#  ```{r, define initial population}
library(ggplot2)
U=0.1 ##mutation rate
N=200000
bottleneck=1/10000 ##dilution
growth=1/bottleneck ## average growth factor
L=3000 ##number of sites

generations=15 ##number of passages
simulations=50

N_pool=20 ##absolute size of the pool to be used for migration simulations (before growth)

#```
#```{r, growth cycles}

Freq_all=(1:generations)*U
Counts_all=c()
N_all=rep(N,generations)

Pool_G=list() ##will contain the genomes of the pooled samples for the migration simulations
Pool_n=list() ##will contain the abundances of the pooled samples for the migration simulations

M=1
for(t in 1:generations){
  Pool_n[[t]]=0
  Pool_G[[t]]=matrix(0,nrow=M,ncol=L)
}

for(pop in 1:simulations){
  print(paste("Population",pop))
  
  M=1 #number of lineages/mutations (SGV)
  id=1:M ## label for each lineage
  
  Genotypes=matrix(0,nrow=M,ncol=L)
  
  frequencies=rep(1/M,M) ##sgv frequencies
  n=N*bottleneck*frequencies ##abundances before growth
  
  
  freq=c() ##will contain the summed frequency of mutations over time
  N_t=c() ##will contain the tot population sizes (after growth) over time
  
  for(t in 1:generations){
    ##growth (expand by x, where x is Poisson distributed with mean = growth)
    N_new=n*rpois(length(n),growth) #vector of new abundances
    N_t=c(N_t,sum(N_new))
    
    #print(sum(N_new))
    ## mutations
    parents=c()
    mut_id=c()
    M_t=M ## temporary number of genotypes
    N_old=N_new
    ## mutation iteration for each genotype
    for(m in 1:M){
      n_i=N_new[m] #abundance of genotype m
      mutations= min(n_i,rpois(1,U*n_i)) # mutations in genotype m
      
      if(mutations>0){
        parents=c(parents,rep(m,mutations)) ##vector of parental id for each mutation
        mut_id=c(mut_id,((M_t+1):(M_t+mutations))) #unique mutation IDs
        ##transform parental cells into new lineages
        N_new[m]=n_i-mutations ## 
        N_new=c(N_new,rep(1,mutations)) ##mutants start with freq 1/N
        M_t=M_t+mutations 
        
      }
    }
    #print(M_t-M)
    id=1:M_t
    
    Counts=c()
    for(j in 1:L){
      Counts=c(Counts,sum(Genotypes[,j]*N_old))
    }
    
    Counts=Counts[which(Counts>0)]
    
    ##compute sum of mutation frequencies: old + new
    freq=c(freq,(sum(Counts)+(M_t-M))/(sum(N_new)))
    
    ##sample genotype for next cycle  (bottleneck)
    N_sample=round(sum(N_new)*bottleneck,digits = 0) ## absolute number of individuals to pick
    sampled_pop=sample(id,N_sample,prob=N_new,replace = T) ## IDs of the sampled genotypes
    id_new=sort(unique(sampled_pop)) #unique IDs
    
    ##compute new abundances before growth
    n=c()
    for(x in id_new){
      n=c(n,length(which(sampled_pop==x)))
    }
    
    ## in order to save memory we update the genome matrix only with those mutations that were still present after the bottleneck
    
    ##update mutant genotypes
    new_mutant=which(id_new>M) ##retrieve mutants positions that survive the sampling
    #mutant_id=id_new[new_mutant]-M ##?
    
    new_mut_id=id_new[which(id_new>M)] ##ID of the mutations
    
    surv=setdiff(id_new,new_mut_id) #survived genotypes from previous cycle
    
    ##update genotype of the de novo mutations
    if(length(new_mut_id)>0){
      G_new=matrix(0,nrow=1,ncol=L) ##will contain the genomes of the new mutants
      for(i in 1:length(new_mut_id)){
        parent=parents[which(mut_id==new_mut_id[i])] ##retrieve parental position in the Genotype matrix
        G=Genotypes[parent,]
        site=sample(1:L,1) # randomly pick a genome site
        G[site]=(1+G[site])%%2 ##mutate it by allowing for back mutations
        G_new=rbind(G_new,G)
      }
    }
    
    ##remove dead genotypes
    if(length(surv)==1){
      Genotypes=t(as.matrix(Genotypes[surv,]))
    }else{Genotypes=Genotypes[surv,]}
    
    
    ##add mutant genomes n the matrix
    if(length(new_mut_id)>0){
      Genotypes=rbind(Genotypes,G_new[2:(length(new_mut_id)+1),])
    }
    
    M=length(n)
    
    ##export a pool to be used in the migration simulations every generations
    id=1:M
    pool=sample(id,N_pool,prob=n,replace = T) 
    id_pool=sort(unique(pool))
    n_pool=c()
    for(i in 1:length(id_pool)){
      n_pool=c(n_pool,length(which(pool==id_pool[i])))
    }
    Pool_n[[t]]=c(Pool_n[[t]],n_pool)
    Pool_G[[t]]=rbind(Pool_G[[t]],Genotypes[id_pool,])
  }
  Counts_all=c(Counts_all,Counts) ##store counts only at last generation
  Freq_all=rbind(Freq_all,freq)
  N_all=rbind(N_all,N_t)
  
}

# Save multiple objects
#save(Counts_all,Freq_all,N_all, file = "Virus_N200000_L3000_B1000_I.RData")
#save(Pool_G,Pool_n, file = "Virus_Pool_N200000_L3000_B1000_I.RData")
