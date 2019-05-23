

library(doSNOW)
cl<-makeCluster(4)
on.exit(stopCluster(cl))
library(foreach)
# do not run the whole function
# or do but it just wont make all the plots just one set of them
Plotmaker <- function(){
  source("InfinitePopFunctions.R")
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  library(viridis)
  
  
  #####run from here #####
  
  # q is for foy, x inv, y inv
  for(q in 1:4){
    # dimensions of plot
    size <- 10
    # dominance on the SAL 1 = male benefit allele dominant
    h1 = 1
    # reduction in recombination
    rfs <-c(.5,1,1,.97)
    rf = rfs[2]
    # the genotype the mutation is introducted into
    intros <- c(10,4,10,10)
    intro <- intros[3]
    # these are the parameters varied on the axis. 
    cd.vals <- seq(.75, 1, length.out = size) # seq(.5,1, length.out = size)
    s.vals <- seq(.1,.6, length.out = size) # seq(0,.5, length.out = size)
    
    # the data stored in this matrix is later plotted. 
    # Ending frequency of the recombination modifying mutation
    result.matrix <- array(dim = c(size,size))
    absmatrix <- array(dim = c(size,size))
    recommatrix <- array(dim = c(size,size))
    premut.eq.geno <- array()
    postmut.eq.geno <- array()
    opts <- list(preschedule=FALSE)
    loop.results <- foreach(i = 1:length(s.vals), .options.multicore=opts,
                            .combine = 'cbind') %doPAR% {
                              
                              result.vector <- c()
                              absmatrix <- array(dim = c(size,size))
                              recommatrix <- array(dim = c(size,size))
                              premut.eq.geno <- array()
                              postmut.eq.geno <- array()
                              for(j in 1:length(s.vals)){
                                
                                cat('\014')
                                cat('abs',q,':',(i - 1) * size + j, 'of' , size^2)
                                cat('\n\n')
                                
                                # going back to reference (no small r)
                                # create a new reference for the beginning of each row
                                if(j == 1 && i == 1){
                                  ##### load some initial genotype frequencies; geno #####
                                  ref.geno <- c(.0625,  # XGARs 
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
                                } else if(j == 1){
                                  ref.geno <- col.ref.geno
                                }
                                geno <- ref.geno
                                
                                # creating parameter vector for this point
                                par  <- c(cd.vals[i], #cd1 
                                          cd.vals[i], #cd2
                                          h1,   #h1
                                          .5,   #h2
                                          0,   #h3
                                          s.vals[j], #sel
                                          .1,  # reduction due to sex geno mismatch
                                          .1,  # r1
                                          .1,  # r2
                                          .1,  # r3
                                          rf,   #rf
                                          .9,  #gs  
                                          .5)   #u   
                                
                                # proceeding to the equilibrium for these params 
                                done <- F
                                while(done == F){
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter=10000,
                                                  prnt = F, 
                                                  tol = 10^-5)[[1]]
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter="stable",
                                                  prnt = F, 
                                                  tol = 10^-5)[[1]]
                                  test.geno <- NewGeno(geno=geno,
                                                       par=par,
                                                       iter=5000,
                                                       prnt = F, 
                                                       tol = 10^-5)[[1]]
                                  if(sum(abs(test.geno - geno)) < 10^(-3)){
                                    done <- T
                                  } 
                                  geno <- test.geno
                                  cat('hi')
                                  rm(test.geno)
                                }
                                rm(done)
                                
                                # saving eq freqs
                                premut.eq.geno[(i - 1) * size + j] <- list(geno)
                                
                                # calculating absolute fitness
                                
                                param <- param.assign(par)
                                pre.mut.abs.fit <- gen(geno, param, return.abs.fit = T)
                                rm(param)
                                
                                # creating new reference for next loop
                                ref.geno <- geno
                                # saving reference for next i loop
                                if(j == 1){
                                  col.ref.geno <- geno
                                }
                                
                                # introduce small r mutation
                                geno[intro] <- .0005
                                
                                # proceeding to equilibrium with the mutation
                                done <- F
                                while(done == F){
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter=10000,
                                                  prnt = F, 
                                                  tol = 10^-3)[[1]]
                                  # proceeding to equilibrium
                                  # start
                                  geno <- NewGeno(geno=geno,
                                                  par=par,
                                                  iter="stable",
                                                  prnt = F, 
                                                  tol = 10^-5)[[1]]
                                  test.geno <- NewGeno(geno=geno,
                                                       par=par,
                                                       iter=5000,
                                                       prnt = F, 
                                                       tol = 10^-5)[[1]]
                                  if(sum(abs(test.geno - geno)) < 10^(-3)){
                                    done <- T
                                  } 
                                  geno <- test.geno
                                  cat('bye')
                                  rm(test.geno)
                                }
                                rm(done)
                                
                                # saving eq freqs
                                postmut.eq.geno[(i - 1) * size + j] <- list(geno)
                                
                                
                                # creating null params for calculating recom load
                                nullpar <- c(cd.vals[i], #cd1 
                                             cd.vals[i], #cd2
                                             h1,   #h1
                                             .5,   #h2
                                             0,   #h3
                                             s.vals[j], #sel
                                             .1,  # reduction due to sex geno mismatch
                                             0,  # r1
                                             0,  # r2
                                             0,  # r3
                                             rf,   #rf
                                             .9,  #gs  
                                             .5)   #u  
                                
                                # these are modified parameter vectors for the gen function
                                nullparam <- param.assign(nullpar)
                                param <- param.assign(par)
                                
                                # this is the freq of the mutation
                                result.matrix[i,j] <- sum(geno[seq(2,32,by = 2)])
                                
                                # this is the recombination load
                                recommatrix[i,j] <- sum(gen(gen(geno,param),param,return.fit = T)) - 
                                  sum(gen(gen(geno,nullparam),nullparam,return.fit = T))
                                
                                # calculating absolute fitness
                                post.mut.abs.fit <- gen(geno, param, return.abs.fit = T)
                                
                                # this is the difference in abs fit before and after the mutation
                                absmatrix[i,j] <- sum(post.mut.abs.fit) - sum(pre.mut.abs.fit)
                                
                                # removing extra values
                                rm(param)
                                rm(nullparam)
                                rm(nullpar)
                                rm(par)
                                rm(geno)
                                
                                
                              }
                            }
    # # if(q==1){
    #   save(result.matrix, file = 'FOYfreq.RData')
    #   save(absmatrix,file = 'FOYabs.RData')
    #   save(recommatrix, file = 'FOYrecom.RData')
    #   save(premut.eq.geno, file = 'FOYpremuteq.RData')
    #   save(postmut.eq.geno, file = 'FOYpostmuteq.RData')
    # } else if(q==2){
    #   save(result.matrix, file = 'XINVfreq.RData')
    #   save(absmatrix,file = 'XINVabs.RData')
    #   save(recommatrix, file = 'XINVrecom.RData')
    #   save(premut.eq.geno, file = 'XINVpremuteq.RData')
    #   save(postmut.eq.geno, file = 'XINVpostmuteq.RData')
    # } else if(q==3){
    #   save(result.matrix, file = 'YINVfreq.RData')
    #   save(absmatrix,file = 'YINVabs.RData')
    #   save(recommatrix, file = 'YINVrecom.RData')
    #   save(premut.eq.geno, file = 'YINVpremuteq.RData')
    #   save(postmut.eq.geno, file = 'YINVpostmuteq.RData')
    # }else if (q==4){
    #   save(absmatrix,file = '97abs.RData')
    #   save(recommatrix, file = '97recom.RData')
    #   save(premut.eq.geno, file = '97premuteq.RData')
    #   save(postmut.eq.geno, file = '97postmuteq.RData')
    # }
  }
  ##### to here #####
  
  
  # ##### Graphing #####
  # load(file = 'FOYfreq.RData')
  # load(file = 'FOYabs.RData')
  # load(file = 'FOYrecom.RData')
  # # this is the frequency plot
  # # creating appropriate data frame
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(result.matrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'freq')
  # new.dat <- as.data.frame(new.dat)
  # 
  # # creating the plot
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = freq, fill = freq) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Frequency of Mutation', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', title = '97% frequency') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(0,0.5,1),labels=c(0,0.5,1),
  #                        limits=c(0,1))
  # 
  # 
  # 
  # # this will be the plot of change in absolute fitness before to after 
  # # the addition of the recombination reducing mutation
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(absmatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'abs')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = abs, fill = abs) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Absolute Fitness', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = '100% rf: Diff in abs fit post - pre mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(absmatrix),max(absmatrix)),
  #                        labels=c(min(absmatrix),max(absmatrix)),
  #                        limits=c(min(absmatrix),max(absmatrix)))
  # 
  # 
  # 
  # # this will be the plot of the recombination load
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(recommatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'recomload')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = recomload, fill = recomload) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Recombination Load', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = '97%: recombination load after mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(recommatrix),0,max(recommatrix)),
  #                        labels=c(min(recommatrix),0,max(recommatrix)),
  #                        limits=c(min(recommatrix),max(recommatrix)))
  # 
  # load(file = 'YINVfreq.RData')
  # load(file = 'YINVabs.RData')
  # load(file = 'YINVrecom.RData')
  # # this is the frequency plot
  # # creating appropriate data frame
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(result.matrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'freq')
  # new.dat <- as.data.frame(new.dat)
  # 
  # # creating the plot
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = freq, fill = freq) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Frequency of Mutation', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', title = 'YINV frequency') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(0,0.5,1),labels=c(0,0.5,1),
  #                        limits=c(0,1))
  # 
  # 
  # 
  # # this will be the plot of change in absolute fitness before to after 
  # # the addition of the recombination reducing mutation
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(absmatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'abs')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = abs, fill = abs) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Absolute Fitness', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = 'YINV: Diff in abs fit post - pre mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(absmatrix),max(absmatrix)),
  #                        labels=c(min(absmatrix),max(absmatrix)),
  #                        limits=c(min(absmatrix),max(absmatrix)))
  # 
  # 
  # 
  # # this will be the plot of the recombination load
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(recommatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'recomload')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = recomload, fill = recomload) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Recombination Load', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = 'YINV: recombination load after mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(recommatrix),0,max(recommatrix)),
  #                        labels=c(min(recommatrix),0,max(recommatrix)),
  #                        limits=c(min(recommatrix),max(recommatrix)))
  # 
  # load(file = 'XINVfreq.RData')
  # load(file = 'XINVabs.RData')
  # load(file = 'XINVrecom.RData')
  # # this is the frequency plot
  # # creating appropriate data frame
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(result.matrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'freq')
  # new.dat <- as.data.frame(new.dat)
  # 
  # # creating the plot
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = freq, fill = freq) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Frequency of Mutation', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', title = 'XINV frequency') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(0,0.5,1),labels=c(0,0.5,1),
  #                        limits=c(0,1))
  # 
  # 
  # 
  # # this will be the plot of change in absolute fitness before to after 
  # # the addition of the recombination reducing mutation
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(absmatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'abs')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = abs, fill = abs) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Absolute Fitness', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = 'XINV: Diff in abs fit post - pre mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(absmatrix),max(absmatrix)),
  #                        labels=c(min(absmatrix),max(absmatrix)),
  #                        limits=c(min(absmatrix),max(absmatrix)))
  # 
  # 
  # 
  # # this will be the plot of the recombination load
  # new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
  #                 rep(seq(.1,.6, length.out = size), each = size),
  #                 as.vector(recommatrix))
  # colnames(new.dat) <- c('corr.det', 'sel.str', 'recomload')
  # new.dat <- as.data.frame(new.dat)
  # 
  # ggplot(new.dat) +
  #   aes(x = sel.str, y = corr.det, z = recomload, fill = recomload) +
  #   geom_tile() + 
  #   coord_equal() +
  #   geom_contour(color = 'white', alpha = .5) +
  #   #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  #   theme_bw() +
  #   labs(fill = 'Recombination Load', x = 'Selection Strength on SAL', 
  #        y = '% Correct Sex Det.', 
  #        title = 'XINV: recombination load after mutation') +
  #   scale_fill_gradientn(colours=heat.colors(100),na.value = "white",
  #                        breaks=c(min(recommatrix),0,max(recommatrix)),
  #                        labels=c(min(recommatrix),0,max(recommatrix)),
  #                        limits=c(min(recommatrix),max(recommatrix)))
  # ##### end graphing#####
  # 
}