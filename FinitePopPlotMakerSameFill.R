library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(viridis)
library(readr)

Ns <- c(300,800,1500, 4000, 5000,10000)
size <- 10
# new.res <- list()
# for(j in 1:15){
#   if(j <= 12){new.res[[j]] <- one[[j]]}
#   if(j > 12){new.res[[j]] <- results[[j-3]]}
# }
mxx <- -100
mnn <- 100
# load('FinitePopResults_N=300_800_1500_4000_5000.RData')
for(j in 1:18){
  if(max(results[[j]]) > mxx){mxx <- max(results[[j]])}
  if(min(results[[j]]) < mnn){mnn <- min(results[[j]])}
}

# this will be the plot of the recombination load


#####
for(i in 1:6){
  N <- Ns[i]
  mx <- mxx - 1/(2 * 300)
  mn <- mnn - 1/(2 * 300)
  # this will be the plot of the recombination load
  new.dat <- results #load('FinitePopResults.RData')
  new.dat <-cbind(rep(seq(10^-8, 10^-3, length.out = size), size),
                  rep(seq(.05, .4, length.out = size), each = size),
                  as.vector(results[[(3 * (i-1) + 1)]] - matrix(1/(2 * N), size, size)))
  colnames(new.dat) <- c('mut.rate', 'rec.dist', 'fix.freq')
  new.dat <- as.data.frame(new.dat)
  
  p <- ggplot(new.dat) +
    aes(x = mut.rate, y =rec.dist, z = fix.freq, fill = fix.freq) +
    geom_tile() + 
    # coord_equal() +
    theme(aspect.ratio=1) +
    geom_contour(color = 'white', alpha = .5) +
    #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
    theme_bw() +
    labs(fill = 'Inversion \nfix freq -\n neutral \nfix freq', x = 'Mutation Rate', 
         y = 'Recombination Distance', 
         title = paste('Freq of Y inv fixation N = ', N)) +
    scale_fill_gradientn(colours=viridis(256),
                         breaks=c(mn,round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7), mx),
                         labels=c(round(mn, digits = 7),
                                  round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7),
                                  round(mx, digits = 7)),
                         limits=c(min(mn,round(min(new.dat$fix.freq), digits = 7)),
                                  max(mx,round(max(new.dat$fix.freq), digits = 7))))
  p
  if(i == 1) a <- p
  if(i == 2) b <- p
  if(i == 3) c <- p
  if(i == 4) h3 <- p
  if(i == 5) h6 <- p
  if(i == 6) h9 <- p
  
  new.dat <-cbind(rep(seq(10^-8, 10^-3, length.out = size), size),
                  rep(seq(.75, 1, length.out = size), each = size),
                  as.vector(results[[(3 * (i-1) + 2)]] - matrix(1/(2 * N), size, size)))
  colnames(new.dat) <- c('mut.rate', 'cd.vals', 'fix.freq')
  new.dat <- as.data.frame(new.dat)
  
  p <- ggplot(new.dat) +
    aes(x = mut.rate, y =cd.vals, z = fix.freq, fill = fix.freq) +
    geom_tile() + 
    # coord_equal() +
    theme(aspect.ratio=1) +
    geom_contour(color = 'white', alpha = .5) +
    #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
    theme_bw() +
    labs(fill = 'Inversion \nfix freq -\n neutral \nfix freq', x = 'Mutation Rate', 
         y = '% Correct Sex Det.', 
         title = paste('Freq of Y inv fixation N = ', N)) +
    scale_fill_gradientn(colours=viridis(256),
                         breaks=c(mn,round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7), mx),
                         labels=c(round(mn, digits = 7),
                                  round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7),
                                  round(mx, digits = 7)),
                         limits=c(min(mn,round(min(new.dat$fix.freq), digits = 7)),
                                  max(mx,round(max(new.dat$fix.freq), digits = 7))))
  p
  if(i == 1) d <- p
  if(i == 2) e <- p
  if(i == 3) f <- p
  if(i == 4) h4 <- p
  if(i ==5) h7 <- p
  if(i == 6) h10 <- p
  
  new.dat <-cbind(rep(seq(.75, 1, length.out = size), size),
                  rep(seq(.05, .4, length.out = size), each = size),
                  as.vector(results[[(3 * (i-1) + 3)]] - matrix(1/(2 * N), size, size)))
  colnames(new.dat) <- c('cd.vals', 'rec.dist', 'fix.freq')
  new.dat <- as.data.frame(new.dat)
  
  p <- ggplot(new.dat) +
    aes(x = cd.vals, y =rec.dist, z = fix.freq, fill = fix.freq) +
    geom_tile() + 
    # coord_equal() +
    theme(aspect.ratio=1) +
    geom_contour(color = 'white', alpha = .5) +
    #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
    theme_bw() +
    labs(fill = 'Inversion \nfix freq -\n neutral \nfix freq', x = '% Correct Sex Determination', 
         y = 'Recombination Distance', 
         title = paste('Freq of Y inv fixation N = ', N)) +
    scale_fill_gradientn(colours=viridis(256),
                         breaks=c(mn,round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7), mx),
                         labels=c(round(mn, digits = 7),
                                  round(min(new.dat$fix.freq), digits = 7),
                                  round(max(new.dat$fix.freq), digits = 7),
                                  round(mx, digits = 7)),
                         limits=c(min(mn,round(min(new.dat$fix.freq), digits = 7)),
                                  max(mx,round(max(new.dat$fix.freq), digits = 7))))
  p
  if(i == 1) g <- p
  if(i == 2) h <- p
  if(i == 3) h2 <- p
  if(i == 4) h5 <- p
  if(i == 5) h8 <- p
  if(i == 6) h11 <- p
}

# a
# b
# c
# h3
# h6
# d
# e
# f
# h4
# h7
# g
# h
# h2
# h5
# h8

grid.arrange(a,
             d,
             g,
             b,
             e,
             h,
             c,
             f,
             h2,
             h3,
             h4,
             h5,
             h6,
             h7,
             h8,
             h9,
             h10,
             h11,
             ncol = 3)
# a
# d
# g
# b
# e
# h
# c
# f
# h2
# h3
# h4
# h5
# h6
# h7
# h8
#####
