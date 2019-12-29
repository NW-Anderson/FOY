# this function performs an entire generation by calling the functions in the correct 
# order mutate, fitness, fitness selection, recombine and reproduce
# it is looped through repeatedly for the forward time simulation
generation <- function(pop,mut.rate,h1,h2,h3,s,t,gs, rectable, cd = .95){
  # We begin each generation by introducing mutations on the general 
  # fitness locus (gfl) on some individuals in the population.
  pop <- mutate(pop, mut.rate, rectable)
  # we then calculate the fitness of every idividual taking into account the gfl
  # the sexually antagonistic locus (sal) and sex
  fits <- c()
  for(i in 1:length(pop)){
    fits <- c(fits, fitness(names(pop)[i],h1,h2,h3,s,t,gs))
  }
  rm(i)
  # the population then undergoes fitness selection. the fit genotypes may increase
  # in number and the unfit may decrease or die off. This function is stochastic 
  # and randomly samples individuals according to fitness
  pop <- fitnessSelection(pop, fits)
  # changing the population from tracking number of individuals to frequency of each genoytype
  geno.freqs <- pop / sum(pop)
  # this next function is where recombination occurs. A certain proportion (corresponding
  # to the probability of such a recombination event occuring) of each genotype is
  # recombined into a new genotype
  haplo.freqs <- ApplyRecTable(geno.freqs, rectable)
  # this function is where new genotypes of thenext generation are created by randomly
  # combining genoypes with probabilities equal to their frequency in the parent population
  # this is another stochastic function
  pop <- Reproduce(rectable, haplo.freqs, sum(pop))
  # returning the next geneeration
  return(pop)
}

# This function takes in a pop and a mutation rate that will turn gfl to g
# the deleterious mutation
# returns the new pop where some individuals have experienced this mutation
mutate <- function(pop, rate, rectable){
  # creating new population
  new.pop <- c(matrix(0, 1, length(pop)))
  names(new.pop) <- names(pop)
  # going through each genotype and mutating the individuals present at the given rate
  for(i in 1:length(pop)){
    # we use twice the rate because there exists two gfl per individual
    # binomial returns how many people of that genotype mutate
    mut.events <- rbinom(1, pop[i], (2 * rate))
    # making sure we are not mutating more individuals than are present
    mut.events <- min(mut.events, pop[i])
    # puts all the individuals who didnt mutate in the new population
    new.pop[i] <- new.pop[i] + pop[i] - mut.events
    # dealing with the mutated individuals
    if(mut.events > 0){
      for(j in 1:mut.events){
        # splitting the genotype into the different haplotypes, sex, and presence 
        # of inversions
        sex <- substr(names(pop)[i],1,1)
        inv <- nchar(names(pop)[i]) == 12
        dubinv <- nchar(names(pop)[i]) == 13
        # because characters are shifted around with inversions
        if(dubinv == T){
          hap1 <- substr(names(pop)[i], 3, 6)
          hap2 <- substr(names(pop)[i], 9, 12)
        }else{
          hap1 <- substr(names(pop)[i], 3, 6)
          hap2 <- substr(names(pop)[i], 8, 11)
        }
        # randomly choose which loci will mutate in this individual and replaces the gfl 
        # with a g (deleterious mutation)
        if(runif(1) < .5){
          hap1 <- paste(substr(hap1,1,1),'g',substr(hap1,3,4), sep = '')
        }else{
          hap2 <- paste(substr(hap2,1,1),'g',substr(hap2,3,4), sep = '')
        }
        # this save is in case the haplotypes are in the incorrect order to match with
        # one of the genotypes present. We cannot simply reverse the order because
        # the inversion may be stuck on the wrong haplotype in the following step
        # if the order is indeed reversed
        hap1save <- hap1
        hap2save <- hap2
        # putting the genotype back together
        # we begin by adding the i to haplotypes if there is an inversion
        if(inv == T){
          hap2 <- paste(hap2,'i', sep = '')
        }else if(dubinv == T){
          hap1 <- paste(hap1,'i', sep = '')
          hap2 <- paste(hap2,'i', sep = '')
        }
        # pasting the sex and haplotypes together to form the genotype
        new.ind.geno <- paste(sex,hap1,hap2)
        # this flips the order of the haplotypes if they are in the wrong order after 
        # the mutation
        if(is.na(match(new.ind.geno, row.names(rectable)))){
          # reversing the order and reverting to the save above so we can stick the 
          # inversions to the correct haplotype
          hap1 <- hap2save
          hap2 <- hap1save
          # adding the inversions to the reversed haplotype order
          if(inv == T){
            hap2 <- paste(hap2,'i', sep = '')
          }else if(dubinv == T){
            hap1 <- paste(hap1,'i', sep = '')
            hap2 <- paste(hap2,'i', sep = '')
          }
          # pasting the sex and haplotypes together to form the genotype
          new.ind.geno <- paste(sex,hap1,hap2)
        }
        ### testing to make sure the reversed order actually exists in the pop table
        if(is.na(match(new.ind.geno,row.names(rectable)))){
          stop('2')
        }
        # putting  the mutated individual in the new genotype
        new.pop[match(new.ind.geno, names(pop))] <- new.pop[match(new.ind.geno, names(pop))] + 1
      }
      ### testing to make sure the nuber of individuals in both the original pop and the newpop
      ### remains constant throughout this process
      if(sum(new.pop) != sum(pop[1:i])){
        stop('3')
      }
    }
  }
  # returning the population where some individuals have experienced the deleterious mutation
  return(new.pop)
}

# takes in a single genotype and selection parameters and returns the fitness
# of that particular genotype
fitness <- function(geno,h1,h2,h3,s,t,gs){
  # breaking apart the genotype into sex, haplotypes and the presence of inversions
  pheno.sex <- substr(geno, 1, 1)
  inv <- nchar(geno) == 12
  dubinv <- nchar(geno) == 13
  # because the char are shifted with the is present in inverted genos
  if(dubinv == T){
    geno.sex <- paste(substr(geno, 3, 3), substr(geno, 9, 9))
    gfl <- paste(substr(geno, 4, 4), substr(geno, 10,10))
    sal <- paste(substr(geno, 5, 5), substr(geno, 11, 11))
    rml <- paste(substr(geno, 6, 6), substr(geno, 12, 12))
  }else{
    geno.sex <- paste(substr(geno, 3, 3), substr(geno, 8, 8))
    gfl <- paste(substr(geno, 4, 4), substr(geno, 9,9))
    sal <- paste(substr(geno, 5, 5), substr(geno, 10, 10))
    rml <- paste(substr(geno, 6, 6), substr(geno, 11, 11))
  }
  # we begin to calculate fitness we begin with fitness =  1 and will mult this by
  # fitness factors corresponding to the fitness of each loci and sex
  fit <- 1
  # we begin by calculating the fitness effects of the sal
  if(sal == 'A A'){
    # the sal has different fitness effects depending on sex
    if(pheno.sex == 'S'){
      fit <- fit * (1+s)
    }else{
      fit <- fit / (1+s)
    }
  }else if(sal == 'A a' | sal == 'a A'){
    if(pheno.sex == 'S'){
      fit <- fit * (1 + h1 * s)
    }else{
      fit <- fit / (1 + h1 * s)
    }
  } # otherwise (sal = a a) the fitness multiplier is 1 so we ignore this case
  # we then calc the fitness effects of the gfl
  # if there is not a deleterious mutation (gfl = G G) the fitness multiplier
  # is one so we ignore this case
  if(gfl == 'G g' | gfl == 'g G'){
    fit <- fit / (1 + h3 * gs)
  }else if(gfl == 'g g'){
    fit <- fit / (1 + gs)
  }
  # we then calculate the fitness effects of any mismatch between the geno 
  # vs pheno sex
  # this is the male without a y
  if(geno.sex == 'X X'){
    if(pheno.sex == 'S'){
      fit <- fit * (1 - t)
    }
  # this is the female with a y
  }else if(geno.sex == 'X Y' | geno.sex == 'Y X'){
    if(pheno.sex == 'D'){
      fit <- fit * (1 - t)
    }
  # this is any individual having too many y
  }else if(geno.sex == 'Y Y'){
    if(pheno.sex == 'S'){
      fit <- fit * (1 - t)
    }else if(pheno.sex == 'D'){
      fit <- fit * (1 - t)^2
    }
  }
  # after multiplyingby the fitness effects we then return the fitness value of this genotype
  return(fit)
}

# this function takes in a population and the fitnesses of every genotype and 
# samples from the population with probabilities proportional to thier fitness
# males and females undergo fitness selection seperately because males compete with 
# males to mate and females compete with other females
### Little questionable if the way the selection is set up is not the culprit
fitnessSelection <- function(pop,fits){
  # seperating the population into males and females
  males <- c()
  females <- c()
  male.fits <- c()
  female.fits <- c()
  for(i in 1 :length(pop)){
    sex <- substr(names(pop)[i],1,1)
    if(sex == 'S'){
      males <- c(males, pop[i])
      male.fits <- c(male.fits, fits[i])
    }else if(sex == 'D'){
      females <- c(females, pop[i])
      female.fits <- c(female.fits, fits[i])
    }
  }
  # sample from the males and females seperately with probabilities proportional to fitness
  # relative to everyone of the same sex
  new.males <- sample(names(males), sum(males), replace = T, prob = (males * male.fits))
  new.females <- sample(names(females), sum(females), replace = T, prob = (females * female.fits))
  # combining males and females into a population
  new.pop <- c(new.males, new.females)
  # turning the list of individulas back into a vector contianing how many individuals exist of each 
  # genotype this vector will not contain any of the genotypes which do not have any individuals present
  # in this population
  new.pop <- c(unlist(table(new.pop)))
  # putting back in the original order and adding the genotypes which have 0 individuals back into the\
  # vector
  reordered.pop <- c()
  for(i in 1:length(pop)){
    if(!is.na(match(names(pop)[i], names(new.pop)))){
      reordered.pop[i] <- new.pop[match(names(pop)[i], names(new.pop))]
    }else{
      reordered.pop[i] <- 0
    }
  }
  names(reordered.pop) <- names(pop)
  # returning the new populaiton
  return(reordered.pop)
}

# takes in a vector of genotype frequencies and returns a vector of haplotype freq
# after recombination events have occured at the given rates
ApplyRecTable <- function(geno.freqs, rectable){
  # this vector will contain the haplotype freqs
  haplo.freqs <- runif(ncol(rectable),0,0)
  names(haplo.freqs) <- colnames(rectable)
  # this will go theough each of the haplotypes and determine what proportion of the population
  # each genotype contributes to that haplotype by multiplying the genotype frequency by the 
  # probability of that genotype recombining to reproduce a given haplotype
  # this is essentially matrix multiplying the vector of genoype freqs by the recomtable 
  for(i in 1:length(geno.freqs)){
    for(j in 1:length(haplo.freqs)){
      # this is essentially summing over each entire column but it in order of rows
      # not by column. its a bit backwards but it is doing the same thing.
      haplo.freqs[j] <- haplo.freqs[j] + geno.freqs[i] * rectable[i,j]
    } 
  }
  # adding haplotype names to the vector of haplotype frequencies
  haplotypes <- c('XGAR', 
                  'XGAr',
                  'XGaR', 
                  'XGar', 
                  'XgAR', 
                  'XgAr',
                  'XgaR', 
                  'Xgar', 
                  'YGAR', 
                  'YGAr', 
                  'YGaR', 
                  'YGar', 
                  'YgAR', 
                  'YgAr', 
                  'YgaR',
                  'Ygar',
                  'YGARi', 
                  'YGAri', 
                  'YGaRi', 
                  'YGari', 
                  'YgARi', 
                  'YgAri', 
                  'YgaRi',
                  'Ygari')
  
  haplotypes <- c(paste('S', haplotypes), paste('D', haplotypes))
  names(haplo.freqs) <- haplotypes
  # returning haplotype frequencies
  return(haplo.freqs)
}

# randomly draws haplotypes correspondng to thier frequencies in the parentpopulation, 
# combines them into individual genotypes and determines sex based on geno sex and 
# probability of sex reversal and then returns the new pop
Reproduce <- function(rectable, haplo.freqs, N, cd = .95){
  # beginning with a vector of 0s which will be filled with the new generation
  pop <- c(matrix(0,1,nrow(rectable)))
  # drawing randomly from the haplotype pool and pasting together
  # drawing haplotypes seperately from sperm and egg
  for(n in 1:N){
    # drawing the first haplotype from sperm corresponding to haplotype freqs in sperm
    hap1 <- sample(names(haplo.freqs)[1:24], 1, 
                   prob = haplo.freqs[1:24])
    # determinging the length of the haplotype this changes depending on whether the 
    # haplotype has been inverted
    crct1 <- nchar(hap1)
    # getting rid of the S or D on the haplotype before we paste the together
    hap1 <- substr(hap1,3,crct1)
    # doing the same with a haplotype drawn from the eggs
    hap2 <- sample(names(haplo.freqs)[25:48], 1, 
                   prob = haplo.freqs[25:48])
    crct2 <- nchar(hap2)
    hap2 <- substr(hap2,3,crct2)
    # pasting the two haplotypes together to form a genotype
    geno <- paste(hap1, 
                  hap2)
    # assigning sex to the new genotype
    # we first determine the genotypic sex
    geno.sex <- paste(substr(geno,1,1),substr(geno,crct1,crct1))
    # we first look at the genotypic sex and then stochastically determine whether the phenotypic 
    # sex will match at a rate corresponding to the correct determination rate (cd)
    if(geno.sex == 'Y Y' | geno.sex == 'X Y' | geno.sex == 'Y X'){
      if(runif(1) < cd){
        geno <- paste('S',geno)
      }else{geno <- paste('D',geno)}
    }else if(geno.sex == 'X X'){
      if(runif(1) < cd){
        geno <- paste('D',geno)
      }else{geno <- paste('S',geno)}
    }
    # if the geno is not one listed it is in the wrong order,
    # reversing the order of the haplotypes in this case
    if(is.na(match(geno, row.names(rectable)))){
      geno <- paste(substr(geno,1,1), hap2, hap1)
    }
    ### testing making sure the changed one matches
    if(is.na(match(geno,row.names(rectable)))){
      stop('4')
    }
    # adding the created individual to the new generation
    pop[match(geno, rownames(rectable))] <- pop[match(geno, rownames(rectable))] + 1
  }
  # adding names to the new pop
  names(pop) <- rownames(rectable)
  # returning the new generation 
  return(pop)
}

# this creates a table with genotypes onhte rows and haplotypes on the column. 
# each ijth entry is the probability of geno i to recombine to produce haplo j
rectablemaker <- function(dst, red.fac, dom.fac, message = F){
  # creating column names (all possible haplotypes)
  haplotypes <- c('XGAR', 
                  'XGAr',
                  'XGaR', 
                  'XGar', 
                  'XgAR', 
                  'XgAr',
                  'XgaR', 
                  'Xgar', 
                  'YGAR', 
                  'YGAr', 
                  'YGaR', 
                  'YGar', 
                  'YgAR', 
                  'YgAr', 
                  'YgaR',
                  'Ygar',
                  'YGARi', 
                  'YGAri', 
                  'YGaRi', 
                  'YGari', 
                  'YgARi', 
                  'YgAri', 
                  'YgaRi',
                  'Ygari')
  
  # creating vector of genotypes (rownames)
  genotypes <- c()
  for(i in 1:24){
    for(j in i:24){
      # every haplotype pairs with itself and every haplotype that 
      # occurs later in the vector (to eliminate reverse repeats such as
      # genotype ij is equivalent to genotype ji where i and j are 
      # haplotypes)
      genotypes <- c(genotypes,paste(haplotypes[i],haplotypes[j]))
    }
  }
  # differentiating between male and female haplotypes
  haplotypes <- c(paste('S', haplotypes), paste('D', haplotypes))
  
  ### dont think this is needed anymore but im scared to get rid of it
  # only dealing with inv on the Y
  # invgeno <- c()
  # doubinvgeno <- c()
  # for(k in 1:length(genotypes)){
  #   if(substr(genotypes[k], 6, 6) == 'Y'){
  #     invgeno <- c(invgeno, paste(genotypes[k],'i', sep = ''))
  #   }
  #   if(substr(genotypes[k], 6, 6) == 'Y' && substr(genotypes[k], 1, 1) == 'Y'){
  #     doubinvgeno <- c(doubinvgeno, paste(paste(substr(genotypes[k], 1, 4), 'i',
  #                                               sep = ''), paste(
  #                                                 substr(genotypes[k], 6, 9), 'i',
  #                                                 sep = '')))
  #   }
  # }
  # genotypes <- c(genotypes, invgeno, doubinvgeno)
  # doubling genotypes and assigning sex
  
  # differentiating between male and female genotypes
  genotypes <- c(paste("S", genotypes), paste("D", genotypes))
  
  # this will contain the recombination probs
  rectable <- array(dim = c(length(genotypes),length(haplotypes)))
  rownames(rectable) <- genotypes
  colnames(rectable) <- haplotypes
  for(k in 1:nrow(rectable)){ # rows
    for(c in 1:ncol(rectable)){ # columns
      target.hap <- substr(haplotypes[c],3,6)
      target.inv <- nchar(haplotypes[c]) == 7
      target.sex <- substr(haplotypes[c],1,1)
      sex <- substr(genotypes[k], 1, 1)
      inv <- nchar(genotypes[k]) == 12
      dubinv <- nchar(genotypes[k]) == 13
      
      hap1 <- substr(genotypes[k], 3, 6)
      hap2 <- substr(genotypes[k], 8, 11)
      if(dubinv == T){
        hap2 <- substr(genotypes[k], 9, 12)
      }
      if(sex == target.sex){
        # there is no inv in the genotype
        if(inv == F && dubinv == F){
          # the target hap is not inverted
          if(target.inv == F){
            # this is the normal case with no inversions.
            # if both of the haplotypes in the genotype are the same...
            if(hap1 == hap2){
              # and if they match the haplotype of interest
              if(hap1 == target.hap){
                # the individual is certain to produce a gamete with the haplo of interest
                val <- 1
                # otherwise it is impossible
              }else{val <- 0}
              # this is the case where interesting things happen with recombination
              # ie the two haplo of a geno do not match
            }else{
              # will contain m m1 m2 mm
              # for match both match the haplo allele, hap1 match hap2 or mismatch 
              # mismatch is when the haplo of interest contains an allele not present in geno
              alvec <- c()
              for(l in 1:4){
                if(substr(target.hap,l,l) == substr(hap1, l, l) &&
                   substr(target.hap,l,l) == substr(hap2, l, l)){
                  alvec <- c(alvec, 'm')
                  # when only the 1st half of the geno has an allele matching the haplo
                }else if(substr(target.hap,l,l) == substr(hap1, l, l)){
                  alvec <- c(alvec, 'm1')
                  # when 2nd geno has an allele matching the haplo of interst
                }else if(substr(target.hap,l,l) == substr(hap2, l, l)){
                  alvec <- c(alvec, 'm2')
                  # when neither do
                }else{alvec <- c(alvec, 'mm')}
              }
              
              
              # FROM HERE DOWN IS PRISTINE - HB
              # calculation recom dst
              if(sex == 'S'){
                if(substr(hap1, 4, 4) == 'R'){
                  if(substr(hap2, 4, 4) == 'R'){
                    # RR
                    r <- dst
                  }else{
                    # Rr
                    r <- dst - dst * dom.fac * red.fac
                  }
                }else{
                  if(substr(hap2, 4, 4) == 'R'){
                    # Rr
                    r <- dst - dst*dom.fac*red.fac
                  }else{
                    # rr
                    r <- dst - dst * red.fac
                  }
                }
              }else if(sex == 'D'){
                r <- dst
              }
              
              # beginning to calculate recombination probs
              # if an allele in the haplo is missing it is impossible to make from geno
              if('mm' %in% alvec){
                val <- 0
                # this is the case where the two halves of the geno share no common alleles
              }else if(!'m' %in% alvec){
                # starting at one half to start multiplying rates because only half the geno is
                # capable of forming the haplo of interest
                val <- .5
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec[i] == alvec[i - 1]){
                    val <- val * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    val <- val * r
                    if(message){cat('r')}
                  }
                }
                # this is the case where two halves of geno have a common allele at a single loci
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 1){
                # if there exists 1 m there are two possible ways to recombine to make
                # the haplotype of interest
                vals <- c(.5,.5)
                # first possible way we will treat the m as an m1
                alvec1 <- alvec
                alvec1[which(alvec %in% 'm')] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec1[i] == alvec1[i - 1]){
                    vals[1] <- vals[1] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[1] <- vals[1] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # the second possible way to recombine we will treat the m as m2
                alvec2 <- alvec
                alvec2[which(alvec %in% 'm')] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec2[i] == alvec2[i - 1]){
                    vals[2] <- vals[2] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[2] <- vals[2] * r
                    if(message){cat('r')}
                  }
                }
                val <- sum(vals)
                # this is the case where the two halves of the geno have two loci in common
                # in this case there are 4 possible ways to recobine and produce the haplo
                # of interest
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 2){
                vals <- c(.5, .5, .5, .5)
                alvec1 <- alvec2 <- alvec3 <- alvec4 <- alvec
                # in the first possibility we treat both m as m1
                alvec1[which(alvec %in% 'm')] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec1[i] == alvec1[i - 1]){
                    vals[1] <- vals[1] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[1] <- vals[1] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # in the second possibility we treat the first m as m1 and the second as m2
                alvec2[which(alvec %in% 'm')[1]] <- 'm1'
                alvec2[which(alvec %in% 'm')[2]] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec2[i] == alvec2[i - 1]){
                    vals[2] <- vals[2] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[2] <- vals[2] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # in the third possibility we treat the first m as m2 and the second as m1
                alvec3[which(alvec %in% 'm')[1]] <- 'm2'
                alvec3[which(alvec %in% 'm')[2]] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec3[i] == alvec3[i - 1]){
                    vals[3] <- vals[3] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[3] <- vals[3] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # in the fourth possibility we both m as m2
                alvec4[which(alvec %in% 'm')] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec4[i] == alvec4[i - 1]){
                    vals[4] <- vals[4] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[4] <- vals[4] * r
                    if(message){cat('r')}
                  }
                }
                val <- sum(vals)
                # this is the case where the two halves of the geno have three loci in common
                # in this case it does not matter if there is a recom or not, one gamete 
                # produced will be the haplotype of interest one will not
                # there is no 4 m case because this was accounted for when the two halves
                # of the geno were equal
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 3){
                val <- .5
              }else(stop(cat('Error',k,c)))
            }
            # the target hap is inverted
          }else if(target.inv == T){
            # this case the geno does not contain an inv but the target hap does
            # impossible
            val <- 0
          }
          # when one of the haplotypes in the geno is inv. No recom should take place
        }else if(inv == T){
          if(target.inv == F){
            # in this case only hap1 can contribute to the target hap because hap2 is inv
            # if there is only one inv
            if(hap1 == target.hap){
              val <- .5
            }else{
              # otherwise hap2 contains the inv and no recom will take place so if hap1
              # is not the same as the target then the target will not be produced
              val <- 0
            }
          }else if(target.inv == T){
            # in this case there is a single inv in the geno so no recom occurs
            # the target haplo has this inv so the only way to produce the target
            # is if hap2 == target
            if(hap2 == target.hap){
              val <- .5
            }else{
              val <- 0
            }
          }
        }else if(dubinv == T){
          # in this case recom will take place in the geno but it will only be able to 
          # produce the target if the target is also inv
          if(target.inv == F){
            val <- 0
          }else{
            # this is same as the case when there are no inv in how it is handled. there
            # exists recom and it will produce haplo with inv as well
            if(hap1 == hap2){
              # and if they match the haplotype of interest
              if(hap1 == substr(target.hap,1,4)){
                # the individual is certain to produce a gamete with the haplo of interest
                val <- 1
                # otherwise it is impossible
              }else{val <- 0}
              # this is the case where interesting things happen with recombination
              # ie the two haplo of a geno do not match
            }else{
              # will contain m m1 m2 mm
              # for match both match the haplo allele, hap1 match hap2 or mismatch 
              # mismatch is when the haplo of interest contains an allele not present in geno
              alvec <- c()
              for(l in 1:4){
                if(substr(target.hap,l,l) == substr(hap1, l, l) &&
                   substr(target.hap,l,l) == substr(hap2, l, l)){
                  alvec <- c(alvec, 'm')
                  # when only the 1st half of the geno has an allele matching the haplo
                }else if(substr(target.hap,l,l) == substr(hap1, l, l)){
                  alvec <- c(alvec, 'm1')
                  # when 2nd geno has an allele matching the haplo of interst
                }else if(substr(target.hap,l,l) == substr(hap2, l, l)){
                  alvec <- c(alvec, 'm2')
                  # when neither do
                }else{alvec <- c(alvec, 'mm')}
              }
              
              
              # FROM HERE DOWN IS PRISTINE - HB
              # calculation recom dst
              if(sex == 'S'){
                if(substr(hap1, 4, 4) == 'R'){
                  if(substr(hap2, 4, 4) == 'R'){
                    # RR
                    r <- dst
                  }else{
                    # Rr
                    r <- dst - dst * dom.fac * red.fac
                  }
                }else{
                  if(substr(hap2, 4, 4) == 'R'){
                    # Rr
                    r <- dst - dst*dom.fac*red.fac
                  }else{
                    # rr
                    r <- dst - dst * red.fac
                  }
                }
              }else if(sex == 'D'){
                r <- dst
              }
              
              # beginning to calculate recombination probs
              # if an allele in the haplo is missing it is impossible to make from geno
              if('mm' %in% alvec){
                val <- 0
                # this is the case where the two halves of the geno share no common alleles
              }else if(!'m' %in% alvec){
                # starting at one half to start multiplying rates because only half the geno is
                # capable of forming the haplo of interest
                val <- .5
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec[i] == alvec[i - 1]){
                    val <- val * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    val <- val * r
                    if(message){cat('r')}
                  }
                }
                # this is the case where two halves of geno have a common allele at a single loci
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 1){
                # if there exists 1 m there are two possible ways to recombine to make
                # the haplotype of interest
                vals <- c(.5,.5)
                # first possible way we will treat the m as an m1
                alvec1 <- alvec
                alvec1[which(alvec %in% 'm')] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec1[i] == alvec1[i - 1]){
                    vals[1] <- vals[1] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[1] <- vals[1] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # the second possible way to recombine we will treat the m as m2
                alvec2 <- alvec
                alvec2[which(alvec %in% 'm')] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec2[i] == alvec2[i - 1]){
                    vals[2] <- vals[2] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[2] <- vals[2] * r
                    if(message){cat('r')}
                  }
                }
                val <- sum(vals)
                # this is the case where the two halves of the geno have two loci in common
                # in this case there are 4 possible ways to recobine and produce the haplo
                # of interest
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 2){
                vals <- c(.5, .5, .5, .5)
                alvec1 <- alvec2 <- alvec3 <- alvec4 <- alvec
                # in the first possibility we treat both m as m1
                alvec1[which(alvec %in% 'm')] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec1[i] == alvec1[i - 1]){
                    vals[1] <- vals[1] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[1] <- vals[1] * r
                    if(message){cat('r')}
                  }
                }
                if(message){cat('\n\n')}
                # in the second possibility we treat the first m as m1 and the second as m2
                alvec2[which(alvec %in% 'm')[1]] <- 'm1'
                alvec2[which(alvec %in% 'm')[2]] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec2[i] == alvec2[i - 1]){
                    vals[2] <- vals[2] * (1 - r)
                    if(message){cat('1-r')}
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[2] <- vals[2] * r
                    if(message){cat('r')}
                  }
                }
                cat('\n\n')
                # in the third possibility we treat the first m as m2 and the second as m1
                alvec3[which(alvec %in% 'm')[1]] <- 'm2'
                alvec3[which(alvec %in% 'm')[2]] <- 'm1'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec3[i] == alvec3[i - 1]){
                    vals[3] <- vals[3] * (1 - r)
                    if(message)cat('1-r')
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[3] <- vals[3] * r
                    if(message)cat('r')
                  }
                }
                cat('\n\n')
                # in the fourth possibility we both m as m2
                alvec4[which(alvec %in% 'm')] <- 'm2'
                for(i in 2:4){
                  # if m1 -> m1 then there cannot be a recom event
                  if(alvec4[i] == alvec4[i - 1]){
                    vals[4] <- vals[4] * (1 - r)
                    if(message)cat('1-r')
                    # if m1 to m2 then recom must occur
                  }else{
                    vals[4] <- vals[4] * r
                    if(message)cat('r')
                  }
                }
                val <- sum(vals)
                # this is the case where the two halves of the geno have three loci in common
                # in this case it does not matter if there is a recom or not, one gamete 
                # produced will be the haplotype of interest one will not
                # there is no 4 m case because this was accounted for when the two halves
                # of the geno were equal
              }else if(table(alvec)[names(table(alvec)) == 'm'] == 3){
                val <- .5
              }else(stop(cat('Error',k,c)))
            }
          }
        }
      }else {val <- 0}
      rectable[k,c] <- val
    }
  }
  return(rectable)
}
