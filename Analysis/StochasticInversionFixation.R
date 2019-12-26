##### parameters  #####
#library(doSNOW)
#cl<-makeCluster(3)
#on.exit(stopCluster(cl))
library(foreach)
library(doMC)
registerDoMC(4)
source('StochasticInternalFunctions.R')
# these will eventually be varied
size <- 10
# noise 
# mu
N.vals <- c(500, 1500, 4000, 10000, 25000)
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
for(n in 1:3){ #length(N.vals)){
  N <- N.vals[n]

  ##### 1st plot #####
  # the first plot is the freq an inversion fixes in the pop plotted against 
  # mutation rate and the recombination distance
  
  # inv fix freqs is the matrix containing inv fix freqs. this is the data for the 1st
  # plot for this pop size
  inv.fix.freqs <- array(dim = c(size, size))
  # we first loop through mutation rates which will be on the x axis
  for(i in 1:length(rates)){
    mut.rate <- rates[i]
    # we then loop through recombination distances between each loci
    for(j in 1:length(rec.dists)){
      # for each new rec dist we must recreate the rectable
      rectable <- rectablemaker(rec.dists[j], 1, .5)
      cat('\014')
      cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
      cat('reaching eq')
      # if this is the first time the rec table is made we go through and find all
      # of the genos without inv and with the FOY recomreducing allele fixed
      if(i == 1 && j ==1){
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
        # to find the equilibrium population we first draw randomly from the orig genos
        # then we run for several hundred generations to allow the pop to reach eq
        # we use large pop numbers so there is less stochasiticity due to sampling
        eq.pop.save <- c() 
        pop <- c(unlist(table(sample(seq(1:length(rownames(rectable)))[orig.genos], 
                                     30000, replace = T))))
        for(x in 1:600){
          if(!is.na(match(x,names(pop)))){
            eq.pop.save[x] <- pop[match(x,names(pop))]
          }else{eq.pop.save[x] <- 0}
        }
        rm(pop)
        names(eq.pop.save) <- rownames(rectable)
        # eq.pop <- sample(rownames(rectable)[orig.genos], 10000, replace = T)
        for(z in 1:250){
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('reaching eq gen:',z)
          eq.pop.save <- generation(eq.pop.save,mut.rate,h1,h2,h3,s,t,gs, rectable)
        }
      }
      eq.pop <- eq.pop.save
      cat('\014')
      cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
      cat('inserting inversions')
      # inv fix will be the count of how many of the k lead to inv fixing in the pop
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
      invfix <- foreach(k = 1:1000, .options.multicore=opts, .combine = 'c') %dopar% {
        # we begin by sampling a finite population from the equil pop we 
        # found before
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
        # done will turn to true when the inv fixes or dies out
        done <- F
        while(done == F){
          # doing 1 generation
          pop <- generation(pop,mut.rate,h1,h2,h3,s,t,gs, rectable)
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          
          # we then calculate the inversion freq in the pop
          # these exist only to find inv freq in each generation
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
          # check if the inv has fixed or died out
          if(inv.freq == 1 | inv.freq == 0){done <- T}
        }
        cat('\n\n')
        # if the inv fixes in the while loop then we add one to the count of times
        # the inversion fixed
        if(inv.freq == 1){TRUE}else{F}
      }
      # the ijth entry will be the proportion of the 1000 tries that the inv fixed
      inv.fix.freqs[i,j] <- sum(invfix)/1000
    }
  }
  # saving the resulting matrix in the results list
  results[[3*(n-1)+1]] <- inv.fix.freqs
  
  
  
  
  ##### 2nd plot ######
  # the second plot plots freq of inv fixation against the mutation rate and the
  # noise in sex determination
  # we begin by creating the rec table for a dist = 1 (it is not varied here)
  rectable <- rectablemaker(.1, 1, .5)
  # we first loop over the rates
  for(i in 1:length(rates)){
    mut.rate <- rates[i]
    # then we loop over the amount of noise in sex determination.
    for(j in 1:length(cd.vals)){
      cd <- cd.vals[j]
      # we first go to eq in a large population 
      # if this is the first time the rec table is made we go through and find all
      # of the genos without inv and with the FOY recomreducing allele fixed
      if(i == 1 && j ==1){
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
        # to find the equilibrium population we first draw randomly from the orig genos
        # then we run for several hundred generations to allow the pop to reach eq
        # we use large pop numbers so there is less stochasiticity due to sampling
        eq.pop.save <- c() 
        pop <- c(unlist(table(sample(seq(1:length(rownames(rectable)))[orig.genos], 
                                     30000, replace = T))))
        for(x in 1:600){
          if(!is.na(match(x,names(pop)))){
            eq.pop.save[x] <- pop[match(x,names(pop))]
          }else{eq.pop.save[x] <- 0}
        }
        rm(pop)
        names(eq.pop.save) <- rownames(rectable)
        # eq.pop <- sample(rownames(rectable)[orig.genos], 10000, replace = T)
        for(z in 1:250){
          cat('\014')
          cat('n = ', n, 'plot 2', 'i =', i ,'j=', j, '\n')
          cat('reaching eq gen:',z)
          eq.pop.save <- generation(eq.pop.save,mut.rate,h1,h2,h3,s,t,gs, rectable)
        }
      }
      eq.pop <- eq.pop.save
      # we then try 1000 times to insert an inv to the pop and see if it fixes
      # or goes extinct
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
      invfix <- foreach(k = 1:1000, .options.multicore=opts, .combine = 'c') %dopar% {
        # first we sample our pop from the eq
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
        # then run to fixation or extinction
        while(done == F){
          pop <- generation(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
          cat('\014')
          cat('n = ', n, 'plot 2', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          
          # these exist only to find inv freq in each generation
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
        # if the inv fixed we add one to the correct count
        if(inv.freq == 1){TRUE}else{F}
      }
      # prop of tries that fixed goes into the ijth element
      inv.fix.freqs[i,j] <- sum(invfix)/1000
    }
  }
  # saving the above matrix into the results list
  results[[3*(n-1)+2]] <- inv.fix.freqs
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
      # if this is the first time the rec table is made we go through and find all
      # of the genos without inv and with the FOY recomreducing allele fixed
      if(i == 1 && j ==1){
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
        # to find the equilibrium population we first draw randomly from the orig genos
        # then we run for several hundred generations to allow the pop to reach eq
        # we use large pop numbers so there is less stochasiticity due to sampling
        eq.pop.save <- c() 
        pop <- c(unlist(table(sample(seq(1:length(rownames(rectable)))[orig.genos], 
                                     50000, replace = T))))
        for(x in 1:600){
          if(!is.na(match(x,names(pop)))){
            eq.pop.save[x] <- pop[match(x,names(pop))]
          }else{eq.pop.save[x] <- 0}
        }
        rm(pop)
        names(eq.pop.save) <- rownames(rectable)
        # eq.pop <- sample(rownames(rectable)[orig.genos], 10000, replace = T)
        for(z in 1:250){
          cat('\014')
          cat('n = ', n, 'plot 3', 'i =', i ,'j=', j, '\n')
          cat('reaching eq gen:',z)
          eq.pop.save <- generation(eq.pop.save,mut.rate,h1,h2,h3,s,t,gs, rectable)
        }
      }
      eq.pop <- eq.pop.save
      # trying 1000 times to fix the inv in the pop
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
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
  