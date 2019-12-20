#### params and population making ####

source('StochasticInternalFunctions.R')
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
N <- 1000 # population size
rate <- .01 # mutation rate
rectable <- rectablemaker(.2, 1, .5)
eq.pop <- c() 
pop <- c(unlist(table(sample(seq(1:length(rownames(rectable))), N, replace = T))))
for(x in 1:600){
  if(!is.na(match(x,names(pop)))){
    eq.pop[x] <- pop[match(x,names(pop))]
  }else{eq.pop[x] <- 0}
}
rm(pop)
rm(x)
names(eq.pop) <- rownames(rectable)

#####

new.pop <- generation(pop = eq.pop, mut.rate = rate, h1,h2,h3,s,t,gs,rectable)

    new.pop <- mutate(pop = eq.pop,rate,rectable) # may have fixed a little somehting or just made it more complicated needlessly but 
    # it chechs out now. left some testing so those may go off if the problem persists
    
    # definitely fixed someting to do with length(mut.events)
    fits <- c()
    for(i in 1:length(eq.pop)){
      fits <- c(fits, fitness(geno = names(new.pop)[i],h1,h2,h3,s,t,gs)) # fixed a parenthesis error with the ts
    }
    rm(i)
    new.pop <- fitnessSelection(pop = new.pop, fits) # this one actually looks okay maybe comb through again
    #creating genotype freqs
    geno.freqs <- new.pop / sum(new.pop)
    
    haplo.freqs <- ApplyRecTable(geno.freqs, rectable) # Im still a little wary of my matrix mult but i think this is okay
    pop <- Reproduce(rectable, haplo.freqs, sum(new.pop))
    