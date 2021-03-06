##### parameters  #####
library(doSNOW)
cl<-makeCluster(3)
on.exit(stopCluster(cl))
library(foreach)
#library(doMC)
#registerDoMC(4)
source('finite.pop.functions.Table.R')
source('finite.pop.generationTable.R')
# these will eventually be varied
size <- 10
# noise 
# mu
N.vals <- c(300, 800, 1500, 4000, 5000, 10000)
# 10-9
rates <- seq(10^-8, 10^-2, length.out = size)
rec.dists <- seq(.05, .4, length.out = size)
cd.vals <- seq(.750, 1, length.out = size) # seq(.5,1, length.out = size)

# dominance factor of Sal
h1 <- 1   #h1
# dom fac of rml
h2 <- .5   #h2
# dom fac of gfl
h3 <- 1   #h3
# selection on sal
s <- .5 #sel
t <- .3 # reduction due to sex geno mismatch
gs <- .5  # selection on gfl

##### parameteers #####

# results will be a list of the matrices of inv freqs that we will later plot
results <- list()
# we create the same 3 plots for different population size. So we first loop through 
# pop size
for(n in 1:length(N.vals)){
  N <- N.vals[n]


        orig.genos <- c()
        for(y in 1:length(rownames(rectable))){
          geno <- rownames(rectable)[y]
          inv <- nchar(geno) == 12
          dubinv <- nchar(geno) == 13
          # because the char are shifted with the is present in inverted genos
          if(dubinv == T){
            rml <- paste(substr(geno, 6, 6), substr(geno, 12, 12))
          }else{
            rml <- paste(substr(geno, 6, 6), substr(geno, 11, 11))
          }
          if(rml == 'r r' && inv == F && dubinv == F){
            orig.genos <- c(orig.genos, y)
          }
        }
  #     }

  ##### 3rd plot #####
  # this plot plots inv fix freq against noise in sex determination and the rec
  # distances between loci
  mut.rate <- rates[5]
  inv.fix.freqs <- array(dim = c(size, size))
  for(i in 1:length(cd.vals)){
    cd <- cd.vals[i]
    for(j in 1:length(rec.dists)){
      # creating the rec table for the appropriate rec dist
      rectable <- rectablemaker(rec.dists[j], 1, .5)
      
      # finding eq pop
      eq.pop <- c() 
      pop <- c(unlist(table(sample(seq(1:length(rownames(rectable)))[orig.genos], 10000, replace = T))))
      for(x in 1:600){
        if(!is.na(match(x,names(pop)))){
          eq.pop[x] <- pop[match(x,names(pop))]
        }else{eq.pop[x] <- 0}
      }
      rm(pop)
      names(eq.pop) <- rownames(rectable)
      
      for(z in 1:250){
        cat('\014')
        cat('n = ', n, 'plot 3', 'i =', i ,'j=', j, '\n')
        cat('reaching eq gen:',z)
        eq.pop <- generation(eq.pop,mut.rate,h1,h2,h3,s,t,gs, rectable)
      }
      # trying 1000 times to fix the inv in the pop
      opts <- list(preschedule = FALSE)
      registerDoSNOW(cl)
      invfix <- foreach(k = 1:1000, .options.multicore=opts, .combine = 'c') %dopar% {
        # sampling our pop from eq 
        pop <- sample(names(eq.pop), N, replace = T, prob = eq.pop) 
        # we then go through and insert an inversion on a random y chromosome in 
        # the pop
        inserted <- F
        for(l in 1:length(pop)){
          geno <- pop[l]
          # because the char are shifted with the is present in inverted genos
          chrom <- substr(geno, 8, 8)
          if(chrom == 'Y' && inserted == F){
            pop[l] <- paste(pop[l], 'i', sep = '')
            inserted <- T
          }
        }
        # we initialize inv.freq to .5 because we exit the while loop when inv.freq
        # reaches either 0 or 1
        inv.freq <- .5
        # count counts the number of generations before we reach fix or extinction
        # only purpose is the print statement
        count <- 1
        
        pop <- c(unlist(table(pop)))
        # putting back in the original order
        reordered.pop <- c()
        for(i in 1:nrow(rectable)){
          if(!is.na(match(rownames(rectable)[i], names(pop)))){
            reordered.pop[i] <- pop[match(rownames(rectable)[i], names(pop))]
          }else{
            reordered.pop[i] <- 0
          }
        }
        pop <- reordered.pop
        names(pop) <- names(eq.pop)
        rm(reordered.pop)
        done <- F
        # running until fixation or extinction
        while(done == F){
          pop <- generation(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
          cat('\014')
          cat('n = ', n, 'plot 3', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          
          # these exist only to find the r freq and inv freq in each generation
          iv <- 0
          y <- 0
          for(m in 1:length(pop)){
            inv <- nchar(names(pop)[m]) == 12
            dubinv <- nchar(names(pop)[m]) == 13
            if(dubinv == T){
              geno.sex <- paste(substr(names(pop)[m], 3, 3), substr(names(pop)[m], 9, 9))
            }else{
              geno.sex <- paste(substr(names(pop)[m], 3, 3), substr(names(pop)[m], 8, 8))
            }
            
            if(geno.sex == 'Y Y'){ y <- y + 2 * pop[m]}
            if(geno.sex == 'Y X' | geno.sex == 'X Y'){y <- y + pop[m]}
            
            if(inv == T){iv <- iv + pop[m]}
            if(dubinv == T){iv <- iv + 2 * pop[m]}
            if(iv > y){stop('wtf', cat(m))}
            
          }
          count <- count + 1
          inv.freq <- iv/y
          if(inv.freq == 1 | inv.freq == 0){done <- T}
        }
        cat('\n\n')
        if(inv.freq == 1){T}else{F}
      }
      inv.fix.freqs[i,j] <- sum(invfix)/1000
    }
  }
  results[[3*(n-1)+3]] <- inv.fix.freqs
}
  