# addapted from recomlaodplotter

RightPlot <- function(){
  source("functions1.R")
  
  #####run from here #####
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  library(foreach)
  library(doMC)
  registerDoMC(2)
  
  for(q in 1:3){
    # dimensions of plot
    size <- 10
    cd = .95
    # reduction in recombination
    rfs <-c(.5,1,1)
    rf = rfs[q]
    # the genotype the mutation is introducted into
    intros <- c(10,4,10)
    intro <- intros[q]
    # these are the parameters varied on the axis. The single point of interest is the single numbers below
    h.vals <- seq(0,1,length.out = size)
    s.vals <- seq(.1,.6, length.out = size) # seq(0,.5, length.out = size)
    
    # the data stored in this matrix is later plotted. Ending frequency of the recombination modifying mutation
    result.matrix <- array(dim = c(size,size))
    devmatrix <- array(dim = c(size,size))
    premut.eq.array <- array()
    postmut.eq.array <- array()
    
    # ***** When you are trying to exame a single point run line by line, the loops will take you through the entire process to get
    # the ending frequencies to be plotted. *****
    
    for(i in 1:length(h.vals)){
      opts <- list(preschedule=FALSE)
      loop.results <- foreach(j = 1:length(s.vals), .options.multicore=opts,
                              .combine = 'cbind') %dopar% {
                                # for(j in 1:length(s.vals)){
                                cat('\014')
                                cat(q,':',(i - 1) * size + j, 'of' , size^2)
                                cat('\n\n')
                                
                                # going back to reference (no small r)
                                # create a new reference for the beginning of each row
                                if(i == 1){
                                  ##### load some initial genotype frequencies; geno #####
                                  ref.geno <- col.ref.geno <- c(.0625,  # XGARs 
                                                                0.000,  # XGArs
                                                                .0625,  # XGaRs
                                                                0.000,  # XGars
                                                                .0625,  # XgARs
                                                                0.000,  # XgArs
                                                                .0625,  # XgaRs
                                                                0.000,  # Xgars
                                                                .0625,  # YGARs
                                                                0.000,  # YGArs
                                                                .0625,  # YGaRs
                                                                0.000,  # YGars
                                                                .0625,  # YgARs
                                                                0.000,  # YgArs
                                                                .0625,  # YgaRs
                                                                0.000,  # Ygars
                                                                #***** EGG *****#
                                                                .125,  # XGARe
                                                                0.00,  # XGAre
                                                                .125,  # XGaRe
                                                                0.00,  # XGare
                                                                .125,  # XgARe
                                                                0.00,  # XgAre
                                                                .125,  # XgaRe
                                                                0.00,  # Xgare
                                                                0.00,  # YGARe
                                                                0.00,  # YGAre
                                                                0.00,  # YGaRe
                                                                0.00,  # YGare
                                                                0.00,  # YgARe
                                                                0.00,  # YgAre
                                                                0.00,  # YgaRe
                                                                0.00)  # Ygare
                                  ##### load some initial genotype frequencies #####
                                } else{
                                  ref.geno <- premut.eq.array[[(i - 2) * size + j]]
                                }
                                if(length(ref.geno) == 32){
                                  geno <- ref.geno
                                } else{
                                    ##### load some initial genotype frequencies; geno #####
                                    geno <- col.ref.geno <- c(.0625,  # XGARs 
                                                                  0.000,  # XGArs
                                                                  .0625,  # XGaRs
                                                                  0.000,  # XGars
                                                                  .0625,  # XgARs
                                                                  0.000,  # XgArs
                                                                  .0625,  # XgaRs
                                                                  0.000,  # Xgars
                                                                  .0625,  # YGARs
                                                                  0.000,  # YGArs
                                                                  .0625,  # YGaRs
                                                                  0.000,  # YGars
                                                                  .0625,  # YgARs
                                                                  0.000,  # YgArs
                                                                  .0625,  # YgaRs
                                                                  0.000,  # Ygars
                                                                  #***** EGG *****#
                                                                  .125,  # XGARe
                                                                  0.00,  # XGAre
                                                                  .125,  # XGaRe
                                                                  0.00,  # XGare
                                                                  .125,  # XgARe
                                                                  0.00,  # XgAre
                                                                  .125,  # XgaRe
                                                                  0.00,  # Xgare
                                                                  0.00,  # YGARe
                                                                  0.00,  # YGAre
                                                                  0.00,  # YGaRe
                                                                  0.00,  # YGare
                                                                  0.00,  # YgARe
                                                                  0.00,  # YgAre
                                                                  0.00,  # YgaRe
                                                                  0.00)  # Ygare
                                    ##### load some initial genotype frequencies #####
                                }
                                
                                par  <- c(cd, #cd1 
                                          cd, #cd2
                                          h.vals[i],   #h1
                                          .5,   #h2
                                          1,   #h3
                                          s.vals[j], #sel
                                          .1,  # reduction due to sex geno mismatch
                                          .1,  # r1
                                          .1,  # r2
                                          .1,  # r3
                                          rf,   #rf
                                          .9,  #gs  #no impact good
                                          .5)   #u   #no impact good
                                
                                # proceeding to the equilibrium for these params 
                                done <- F
                                while(done == F){
                                  if(length(geno) != 32){stop('1')}
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter=10000,
                                                  prnt = F, 
                                                  tol = 10^-5)[[1]]
                                  # #geno <- NewGeno(geno=geno,
                                  #                 par=par,
                                  #                 iter="stable",
                                  #                 prnt = F, 
                                  #                 tol = 10^-5)[[1]]
                                  if(length(geno) != 32){stop('2')}
                                  test.geno <- NewGeno(geno=geno,
                                                       par=par,
                                                       iter=5000,
                                                       prnt = F, 
                                                       tol = 10^-5)[[1]]
                                  if(sum(abs(test.geno - geno)) < 10^(-3)){
                                    done <- T
                                  } 
                                  done <- T
                                  geno <- test.geno
                                  cat('hi')
                                  rm(test.geno)
                                }
                                rm(done)
                                
                                # saving eq freqs
                                premut.eq.geno <- list(geno)
                                
                                # creating new reference
                                ref.geno <- geno
                                if(j == 1){
                                  col.ref.geno <- geno
                                }
                                
                                # introduce small r mutation
                                geno[intro] <- .0005
                                
                                # proceeding to the equilibrium for these params 
                                done <- F
                                while(done == F){
                                  if(length(geno) != 32){stop('3')}
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter=10000,
                                                  prnt = F, 
                                                  tol = 10^-5)[[1]]
                                  # geno <- NewGeno(geno=geno,
                                  #                 par=par,
                                  #                 iter="stable",
                                  #                 prnt = F, 
                                  #                 tol = 10^-5)[[1]]
                                  if(length(geno) != 32){stop('4')}
                                  test.geno <- NewGeno(geno=geno,
                                                       par=par,
                                                       iter=5000,
                                                       prnt = F, 
                                                       tol = 10^-5)[[1]]
                                  if(sum(abs(test.geno - geno)) < 10^(-3)){
                                    done <- T
                                  } 
                                  done <- T
                                  geno <- test.geno
                                  cat('bye')
                                  rm(test.geno)
                                }
                                rm(done)
                                
                                # saving eq freqs
                                postmut.eq.geno <- list(geno)
                                
                                nullpar  <- c(1, #cd1 
                                              1, #cd2
                                              h.vals[i],   #h1
                                              .5,   #h2
                                              1,   #h3
                                              s.vals[j], #sel
                                              .1,  # reduction due to sex geno mismatch
                                              .1,  # r1
                                              .1,  # r2
                                              .1,  # r3
                                              rf,   #rf
                                              .9,  #gs  #no impact good
                                              .5)   #u   #no impact good
                                
                                param <- param.assign(par)
                                nullparam <- param.assign(nullpar)
                                
                                results <- list(sum(geno[seq(2,32,by = 2)]), sum(gen(gen(geno,param), param, return.fit = T)) - 
                                                  sum(gen(gen(geno,nullparam), nullparam, return.fit = T)), premut.eq.geno,
                                                postmut.eq.geno)
                                names(results) <- c('freq', 'detload', 'preeq', 'posteq')
                                results
                              }
      for(k in 1:(length(loop.results)/4)){
        result.matrix[i,k] <- as.numeric(loop.results[[1 + 4*(k - 1)]])
        devmatrix[i,k] <- as.numeric(loop.results[[2 + 4*(k - 1)]])
        premut.eq.array <- cbind(premut.eq.array, loop.results[[3 +  4*(k-1)]])
        postmut.eq.array <- cbind(postmut.eq.array, loop.results[[4 +  4*(k-1)]])
      }
      cat('row:', i)
    }
    if(q==1){
      save(result.matrix, file = 'FOYfreq.RData')
      save(devmatrix,file = 'FOYdev.RData')
      save(premut.eq.geno, file = 'FOYpremuteq.RData')
      save(postmut.eq.geno, file = 'FOYpostmuteq.RData')
    } else if(q==2){
      save(result.matrix, file = 'XINVfreq.RData')
      save(devmatrix,file = 'XINVdev.RData')
      save(premut.eq.geno, file = 'XINVpremuteq.RData')
      save(postmut.eq.geno, file = 'XINVpostmuteq.RData')
    } else if(q==3){
      save(result.matrix, file = 'YINVfreq.RData')
      save(devmatrix,file = 'YINVdev.RData')
      save(premut.eq.geno, file = 'YINVpremuteq.RData')
      save(postmut.eq.geno, file = 'YINVpostmuteq.RData')
    }
  }
  ##### to here #####
  
  ##### Graphing #####
  # this is the frequency plot
  # creating appropriate data frame
  new.dat <-cbind(rep(seq(0,1, length.out = size), size), 
                  rep(seq(.05,.6, length.out = size), each = size),
                  as.vector(result.matrix))
  colnames(new.dat) <- c('dom.fac', 'sel.str', 'freq')
  new.dat <- as.data.frame(new.dat)
  
  # creating the plot
  ggplot(new.dat) +
    aes(x = sel.str, y = dom.fac, z = freq, fill = freq) +
    geom_tile() + 
    coord_equal() +
    geom_contour(color = 'white', alpha = .5) +
    #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
    theme_bw() +
    labs(fill = 'Frequency of Mutation', x = 'Selection Strength on SAL', 
         y = 'SAL Dominance Factor', title = 'FOY') +
    scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
                         breaks=c(min(result.matrix), max(result.matrix)),labels=c(min(result.matrix), max(result.matrix)),
                         limits=c(min(result.matrix), max(result.matrix)))
  #save(result.matrix,file = 'X.ns.Arec.radd.rf1.RData')
  
  # this will be the plot of change in recom load from before to after the addition of the recombination reducing mutation
  new.dat <-cbind(rep(seq(0,1, length.out = size), size), 
                  rep(seq(.05,.6, length.out = size), each = size),
                  as.vector(devmatrix))
  colnames(new.dat) <- c('dom.fac', 'sel.str', 'dev')
  new.dat <- as.data.frame(new.dat)
  
  ggplot(new.dat) +
    aes(x = sel.str, y = dom.fac, z = dev, fill = dev) +
    geom_tile() + 
    coord_equal() +
    geom_contour(color = 'white', alpha = .5) +
    #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
    theme_bw() +
    labs(fill = 'Diff in fit', x = 'Selection Strength on SAL', 
         y = 'SAL Dominance Factor', title = 'FOY: Diff in fit with and without noise') +
    scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
                         breaks=c(min(devmatrix),max(devmatrix)),labels=c(min(devmatrix),max(devmatrix)),
                         limits=c(min(devmatrix),max(devmatrix)))

  ##### end graphing#####
  
}