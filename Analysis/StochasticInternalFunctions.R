# This function takes in a pop and a rate and at the given rate it will turn gfl to g 
# returns the new pop
mutate <- function(pop, rate, rectable){
  # creating new population
  new.pop <- c(matrix(0, 1, length(pop)))
  names(new.pop) <- names(pop)
  # going through each genotype type
  for(i in 1:length(pop)){
    # we use twice the rate because there exists two gfl per individual
    # binomial returns how many people of that genotype mutate
    mut.events <- rbinom(1, pop[i], (2 * rate))
    mut.events <- min(pop[i], mut.events)
    # puts all the individuals who didnt mutate in the new population
    new.pop[i] <- new.pop[i] + pop[i] - mut.events
    # dealing with the mutated individuals
    if(mut.events > 0){
      sex <- substr(names(pop)[i],1,1)
      inv <- nchar(names(pop)[i]) == 12
      dubinv <- nchar(names(pop)[i]) == 13
      # because characters are shifted around with inversions
      if(dubinv == T){
        hap1 <- substr(names(pop)[i], 3, 6)
        hap2 <- substr(names(pop)[i], 9, 12)
      }else if(inv == T){
        hap1 <- substr(names(pop)[i], 3, 6)
        hap2 <- substr(names(pop)[i], 8, 11)
      }else{
        hap1 <- substr(names(pop)[i], 3, 6)
        hap2 <- substr(names(pop)[i], 8, 11)
      }
      for(j in 1:length(mut.events)){
        # randomly choose which loci will mutate in this individual
        if(runif(1) < .5){
          hap1 <- paste(substr(hap1,1,1),'g',substr(hap1,3,4), sep = '')
        }else{
          hap2 <- paste(substr(hap2,1,1),'g',substr(hap2,3,4), sep = '')
        }
        # putting the genotype back together
        if(inv == T){
          hap2 <- paste(hap2,'i', sep = '')
        }else if(dubinv == T){
          hap1 <- paste(hap1,'i', sep = '')
          hap2 <- paste(hap2,'i', sep = '')
        }
        new.ind.geno <- paste(sex,hap1,hap2)
        # this flips the order of the haplotypes if they are in the wrong order after 
        # the mutation
        if(is.na(match(new.ind.geno, row.names(rectable)))){
          new.ind.geno <- paste(sex,hap2,hap1)
        }
        # putting  the mutated individual in the new genotype
        new.pop[match(new.ind.geno, names(pop))] <- new.pop[match(new.ind.geno, names(pop))] + 1
      }
    }
  }
  return(new.pop)
}



# takes in a single genotype and selection parameters and returns the fitness
# of the genotype
fitness <- function(geno,h1,h2,h3,s,t,gs){
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
  # beginning to calculate fitness we begin by 1 and will mult values by it
  fit <- 1
  # we begin by calculating the fitness effects of the sal
  if(sal == 'A A'){
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
  } # otherwise the fitness effect is 1 so we ignore
  # we then calc the fitness effects of the gfl
  if(gfl == 'G g' | gfl == 'g G'){
    fit <- fit / (1 + h3 * gs)
  }else if(gfl == 'g g'){
    fit <- fit / (1 + gs)
  }
  # then the fitness effects of the geno vs pheno sex
  if(geno.sex == 'X X'){
    if(pheno.sex == 'S'){
      fit <- fit * 1 - t
    }
  }else if(geno.sex == 'X Y' | geno.sex == 'Y X'){
    if(pheno.sex == 'D'){
      fit <- fit * 1 - t
    }
  }else if(geno.sex == 'Y Y'){
    if(pheno.sex == 'S'){
      fit <- fit * 1 - t
    }else if(pheno.sex == 'D'){
      fit <- fit * (1 - t)^2
    }
  }
  return(fit)
}

# this function takes in a population and a vector of fitnesses and samples from the
# population with probabilities proportional to thier relative fitness
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
  # sample from the males and females with probabilities proportional to fitness
  new.males <- sample(names(males), sum(males), replace = T, prob = (males * male.fits))
  new.females <- sample(names(females), sum(females), replace = T, prob = (females * female.fits))
  # combining to return
  new.pop <- c(new.males, new.females)
  new.pop <- c(unlist(table(new.pop)))
  # putting back in the original order
  reordered.pop <- c()
  for(i in 1:length(pop)){
    if(!is.na(match(names(pop)[i], names(new.pop)))){
    reordered.pop[i] <- new.pop[match(names(pop)[i], names(new.pop))]
    }else{
      reordered.pop[i] <- 0
    }
  }
  names(reordered.pop) <- names(pop)
  return(reordered.pop)
}

# takes in a vector of genotype frequencies and returns a vector of haplotype freq
# based on recombination rate of the supplied genotypes
ApplyRecTable <- function(geno.freqs, rectable){
  # this vector will contain the haplotype freq
  haplo.freqs <- runif(ncol(rectable),0,0)
  names(haplo.freqs) <- colnames(rectable)
  for(i in 1:length(geno.freqs)){
    for(j in 1:length(haplo.freqs)){
      # this is essentially summing over an entire column
      haplo.freqs[j] <- haplo.freqs[j] + geno.freqs[i] * rectable[i,j]
    } 
  }
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
  return(haplo.freqs)
}

# draws haplotypes at the given frequencies, combines them and calculates sex
# and returns a vector of genotypes that make up the individuals in a pop
Reproduce <- function(rectable, haplo.freqs, N, cd = .95){
  pop <- c(matrix(0,1,nrow(rectable)))
  # drawing randomly from the haplotype pool and pasting together
  # drawing seperately from haplotype from sperm and egg
  for(n in 1:N){
    hap1 <- sample(names(haplo.freqs)[1:24], 1, 
                          prob = haplo.freqs[1:24])
    crct1 <- nchar(hap1)
    hap1 <- substr(hap1,3,crct1)
    hap2 <- sample(names(haplo.freqs)[25:48], 1, 
                          prob = haplo.freqs[25:48])
    crct2 <- nchar(hap2)
    hap2 <- substr(hap2,3,crct2)
    geno <- paste(hap1, 
                  hap2)
    # assigning sex to the new genotype
    geno.sex <- paste(substr(geno,1,1),substr(geno,crct1,crct1))
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
    pop[match(geno, rownames(rectable))] <- pop[match(geno, rownames(rectable))] + 1
  }
  names(pop) <- rownames(rectable)
  return(pop)
}

# this creates a table with genotypes onhte rows and haplotypes on the column. 
# each ijth entry is the probability of geno i to recombine to produce haplo j
rectablemaker <- function(dst, red.fac, dom.fac, message = F){
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
  
  # creating vector of genotypes
  genotypes <- c()
  for(i in 1:24){
    for(j in i:24){
      # every haplotype pairs with itself and every haplotype that 
      # occurs later in the vector to eliminate reverse repeats
      genotypes <- c(genotypes,paste(haplotypes[i],haplotypes[j]))
    }
  }
  haplotypes <- c(paste('S', haplotypes), paste('D', haplotypes))
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
  genotypes <- c(paste("S", genotypes), paste("D", genotypes))
  
  # this will contain the recombination probs
  rectable <- array(dim = c(length(genotypes),length(haplotypes)))
  rownames(rectable) <- genotypes
  colnames(rectable) <- haplotypes
  # what recombination values are we using here? r1-9? rec1-3? diff for sexes?
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
                  if(substr(hap2, 2, 2) == 'R'){
                    # RR
                    r <- dst
                  }else{
                    # Rr
                    r <- dst - dst * dom.fac * red.fac
                  }
                }else{
                  if(substr(hap2, 2, 2) == 'R'){
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
                  if(substr(hap2, 2, 2) == 'R'){
                    # RR
                    r <- dst
                  }else{
                    # Rr
                    r <- dst - dst * dom.fac * red.fac
                  }
                }else{
                  if(substr(hap2, 2, 2) == 'R'){
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