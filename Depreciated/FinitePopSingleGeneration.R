
# this function performs an entire generation by calling the functions in the correct 
# order mutate, fitness, fitness selection, recombine and reproduce
generation <- function(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95){
  # we first mutate the gfl on some ind in the pop
  pop <- mutate(pop, mut.rate, rectable)

  for(i in 1:length(pop)){
    if(is.na(match(pop[i], row.names(rectable)))){stop(cat('NA in for(i in 1:n)',i))}
  }
  fits <- c()
  for(i in 1:length(pop)){
    fits <- c(fits, fitness(pop[i],h1,h2,h3,s,t,gs))
  }
  rm(i)
  pop <- fitnessSelection(pop, fits)
  pop <- c(pop[[1]], pop[[2]])
  #creating genotype freqs
  geno.freqs <- runif(nrow(rectable),0,0)
  names(geno.freqs) <- row.names(rectable)
  for(i in 1:length(pop)){
    # this step adds a single individual the geno freqs each time it occurs in pop
    geno.freqs[match(pop[i], row.names(rectable))] <- 
      geno.freqs[match(pop[i], row.names(rectable))] + 1
    if(is.na(geno.freqs[match(pop[i], row.names(rectable))])){
      stop(cat('NA in for(i in 1:n)',i))
    }
    if(is.na(match(pop[i], row.names(rectable)))){stop(cat('NA in for(i in 1:n)',i))}
  }
  # dividing by pop size to turn into freqs
  geno.freqs <- geno.freqs / length(pop)
  haplo.freqs <- ApplyRecTable(geno.freqs, rectable)
  pop <- Reproduce(rectable, haplo.freqs, length(pop))
  return(pop)
}