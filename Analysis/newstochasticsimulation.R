##### parameters  #####
#library(doSNOW)
#cl<-makeCluster(3)
#on.exit(stopCluster(cl))
library(foreach)
library(doMC)
registerDoMC(2)
source('StochasticInternalFunctions.R')
size <- 3
# noise 
# mu
N.vals <- c(800, 4000, 10000, 25000)
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
for(n in 1:4){ 
  N <- N.vals[n]
  
  ##### INVERSION FIX FREQ ~ MUTATION RATE + RECOMBINATION DISTANCE #####
  # the first plot is the freq an inversion fixes in the pop plotted against 
  # mutation rate and the recombination distance
  
  # inv fix freqs is the matrix containing inv fix freqs. this is the data for the 1st
  # plot for this pop size
  inv.fix.freqs <- array(dim = c(size, size,200))
  # we loop through every square on the graph we do every column inside every row so we loop through
  # the rows (rec dist) first
  for(i in 1:size){
    # for each new rec dist we must recreate the rectable
    rectable <- rectablemaker(rec.dists[i], 1, .5)
    # we then loop through every column (mut rate) within each loop of the rows
    for(j in 1:size){
      mut.rate <- rates[j]
      # creating print statement
      cat('\014')
      cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
      cat('reaching eq')
      # if this is the first time the rec table is made we go through and find all
      # of the genos without inv and with the FOY recomreducing allele fixed
      if(i == 1 && j ==1){
        orig.genos <- DetermineOrigGenos(rectable)
        # to find the equilibrium population we first draw randomly from the orig genos
        # then we run for several hundred generations to allow the pop to reach eq
        # we use large pop numbers so there is less stochasiticity due to sampling
        eq.pop.save <- MakeRandomPop(rectable, orig.genos, N = 30000)
        # running the random pop to equilibrium
        for(z in 1:50){ #0){
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('reaching eq save gen:',z)
          eq.pop.save <- generation2.0(eq.pop.save,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95)
        }
      }
      # to find the equilibrium population we first draw randomly from the orig genos
      # then we run for several hundred generations to allow the pop to reach eq
      # we use large pop numbers so there is less stochasiticity due to sampling
      eq.pop <- eq.pop.save
      for(z in 1:25){ #0){
        cat('\014')
        cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
        cat('reaching eq gen:',z)
        eq.pop <- generation2.0(eq.pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95)
      }
      # before we go and insert inversions we are going to reduce our oppulation down to the correct size 
      eq.pop <- SampleEqPop(rectable, eq.pop, N)
      cat('\014')
      cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
      cat('inserting inversions')
      # inv fix will be the count of how many of the k lead to inv fixing in the pop
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
      invfix <- foreach(k = 1:200, .options.multicore=opts, .combine = 'c') %dopar% { #1000, .options.multicore=opts, .combine = 'c') %dopar% {
        pop <- eq.pop
        # we then go through and insert an inversion on a random y chromosome in 
        # the pop
        pop <- InsertInversion(pop)
        # after we have inserted the inversion our next objective is to simulate forward 
        # in time until the inversion either fixes or dies out in the population
        # we initialize inv.freq to .5 because we exit the while loop when inv.freq
        # reaches either 0 or 1
        inv.freq <- .5
        # count counts the number of generations before we reach fix or extinction
        # only purpose is the print statement
        count <- 1
        # done will turn to true when the inv fixes or dies out. it controls the while loop
        done <- F
        while(done == F){
          # each loop is doing one generation
          pop <- generation2.0(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95)
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          # every generation we calculate the inversion freq in the pop
          inv.freq <- CalcInvFreq(pop)
          count <- count + 1
          # check if the inv has fixed or died out
          if(inv.freq == 1 | inv.freq == 0 | count > (4 * N)){done <- T}
        }
        cat('\n\n')
        # if the inv fixes in the while loop then we add one to the count of times
        # the inversion fixed
        inv.freq
      }
      # the ijth entry will be the proportion of the 1000 tries that the inv fixed
      inv.fix.freqs[i,j,] <- invfix
    }
  }
  # saving the resulting matrix in the results list
  results[[3*(n-1)+1]] <- inv.fix.freqs
  
  
  
  
  ##### INVERSION FIX FREQ ~ MUTATION RATE + NOISE IN SEX DETERMINATION ######
  # the second plot plots freq of inv fixation against the mutation rate and the
  # noise in sex determination
  # we begin by creating the rec table for a dist = 1 (it is not varied here)
  rectable <- rectablemaker(.1, 1, .5)
  # we loop through rows (the percent correct sexual development) first
  for(i in 1:size){
    cd <- cd.vals[i]
    # we then loop through the columns or the mutatin rate
    for(j in 1:size){
      mut.rate <- rates[j]
      # we then go to equilibrium for a large population from the same saved
      # population from earlier this keeps the plots from having wildly 
      # different starting populations
      eq.pop <- eq.pop.save
      for(z in 1:250){
        cat('\014')
        cat('n = ', n, 'plot 2', 'i =', i ,'j=', j, '\n')
        cat('reaching eq gen:',z)
        eq.pop <- generation2.0(eq.pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
      }
      eq.pop <- SampleEqPop(rectable, eq.pop, N)
      cat('\014')
      cat('n = ', n, 'plot 2', 'i =', i ,'j=', j, '\n')
      cat('inserting inversions')
      # inv fix will be the count of how many of the k lead to inv fixing in the pop
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
      invfix <- foreach(k = 1:200, .options.multicore=opts, .combine = 'c') %dopar% {
        pop <- eq.pop
        # we then go through and insert an inversion on a random y chromosome in 
        # the pop
        pop <- InsertInversion(pop)
        # after we have inserted the inversion our next objective is to simulate forward 
        # in time until the inversion either fixes or dies out in the population
        # we initialize inv.freq to .5 because we exit the while loop when inv.freq
        # reaches either 0 or 1
        inv.freq <- .5
        # count counts the number of generations before we reach fix or extinction
        # only purpose is the print statement
        count <- 1
        # done will turn to true when the inv fixes or dies out. it controls the while loop
        done <- F
        while(done == F){
          # each loop is doing one generation
          pop <- generation2.0(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          # every generation we calculate the inversion freq in the pop
          inv.freq <- CalcInvFreq(pop)
          count <- count + 1
          # check if the inv has fixed or died out
          if(inv.freq == 1 | inv.freq == 0 | count > (4 * N)){done <- T}
        }
        cat('\n\n')
        # if the inv fixes in the while loop then we add one to the count of times
        # the inversion fixed
        inv.freq
      }
      # the ijth entry will be the proportion of the 1000 tries that the inv fixed
      inv.fix.freqs[i,j,] <- invfix
    }
  }
  # saving the above matrix into the results list
  results[[3*(n-1)+2]] <- inv.fix.freqs
  
  
  
  
  ##### INVERSION FIX FREQ ~ NOISE IN SEX DETERMINATION + RECOMBINATION DISTANCE #####
  # this plot plots inv fix freq against noise in sex determination and the rec
  # distances between loci
  mut.rate <- rates[2]
  inv.fix.freqs <- array(dim = c(size, size))
  # we first loop through the rows (recombination distance)
  for(i in 1:size){
    # creating the rec table for the appropriate rec dist
    rectable <- rectablemaker(rec.dists[i], 1, .5)
    # we then loop through the columns (percent correct sex determinatin)
    for(j in 1:size){
      cd <- cd.vals[j]
      # we then go to equilibrium for a large population from the same saved
      # population from earlier this keeps the plots from having wildly 
      # different starting populations
      eq.pop <- eq.pop.save
      for(z in 1:250){
        cat('\014')
        cat('n = ', n, 'plot 3', 'i =', i ,'j=', j, '\n')
        cat('reaching eq gen:',z)
        eq.pop <- generation2.0(eq.pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
      }
      eq.pop <- SampleEqPop(rectable, eq.pop, N)
      cat('\014')
      cat('n = ', n, 'plot 3', 'i =', i ,'j=', j, '\n')
      cat('inserting inversions')
      # inv fix will be the count of how many of the k lead to inv fixing in the pop
      opts <- list(preschedule = FALSE)
      # registerDoSNOW(cl)
      invfix <- foreach(k = 1:1000, .options.multicore=opts, .combine = 'c') %dopar% {
        pop <- eq.pop
        # we then go through and insert an inversion on a random y chromosome in 
        # the pop
        pop <- InsertInversion(pop)
        # after we have inserted the inversion our next objective is to simulate forward 
        # in time until the inversion either fixes or dies out in the population
        # we initialize inv.freq to .5 because we exit the while loop when inv.freq
        # reaches either 0 or 1
        inv.freq <- .5
        # count counts the number of generations before we reach fix or extinction
        # only purpose is the print statement
        count <- 1
        # done will turn to true when the inv fixes or dies out. it controls the while loop
        done <- F
        while(done == F){
          # each loop is doing one generation
          pop <- generation2.0(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd)
          cat('\014')
          cat('n = ', n, 'plot 1', 'i =', i ,'j=', j, '\n')
          cat('k = ', k, 'gen' = count)
          # every generation we calculate the inversion freq in the pop
          inv.freq <- CalcInvFreq(pop)
          count <- count + 1
          # check if the inv has fixed or died out
          if(inv.freq == 1 | inv.freq == 0 | count > (4 * N)){done <- T}
        }
        cat('\n\n')
        # if the inv fixes in the while loop then we add one to the count of times
        # the inversion fixed
        inv.freq
      }
      inv.fix.freqs[i,j,] <- invfix
    }
  }
  # saving the above matrix into the results list
  results[[3*(n-1)+3]] <- inv.fix.freqs
}