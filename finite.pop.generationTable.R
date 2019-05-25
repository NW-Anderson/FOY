
# this function performs an entire generation by calling the functions in the correct 
# order mutate, fitness, fitness selection, recombine and reproduce
generation <- function(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95){
  # we first mutate the gfl on some ind in the pop
  pop <- mutate(pop, mut.rate, rectable)
  # calculating fitnesses
  fits <- c()
  for(i in 1:length(pop)){
    fits <- c(fits, fitness(names(pop)[i],h1,h2,h3,s,t,gs))
  }
  rm(i)
  pop <- fitnessSelection(pop, fits)
  #creating genotype freqs
  geno.freqs <- pop / sum(pop)
  
  haplo.freqs <- ApplyRecTable(geno.freqs, rectable)
  pop <- Reproduce(rectable, haplo.freqs, sum(pop))
  return(pop)
}