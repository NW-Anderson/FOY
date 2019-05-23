library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(viridis)
library(readr)

Ns <- c(300,800,1500, 4000)
size <- 10
i = 4

mx <- -100
mn <- 100
# load('FinitePopResults.RData')
for(j in 1:12){
  if(max(results[[j]]) > mx){mx <- max(results[[j]])}
  if(min(results[[j]]) < mn){mn <- min(results[[j]])}
}
N <- Ns[i]
mx <- mx - 1/(2 * N)
mn <- mn - 1/(2 * N)
# this will be the plot of the recombination load


#####
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
                       breaks=c(mn,round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5), mx),
                       labels=c(round(mn, digits = 7),
                                round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5),
                                round(mx, digits = 7)),
                       limits=c(mn,mx))
p


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
                       breaks=c(mn,round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5), mx),
                       labels=c(round(mn, digits = 7),
                                round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5),
                                round(mx, digits = 7)),
                       limits=c(mn,mx))
p



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
                       breaks=c(mn,round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5), mx),
                       labels=c(round(mn, digits = 7),
                                round(min(new.dat$fix.freq), digits = 5),
                                round(max(new.dat$fix.freq), digits = 5),
                                round(mx, digits = 7)),
                       limits=c(mn,mx))
p
#####
