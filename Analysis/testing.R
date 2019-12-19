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
names(eq.pop) <- rownames(rectable)
