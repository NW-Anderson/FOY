library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(viridis)
library(readr)


size <- 10

# this will be the plot of the recombination load

# new.dat <-cbind(rep(seq(.75,1, length.out = size), size), 
#                 rep(seq(.1,.6, length.out = size), each = size),
#                 as.vector(recommatrix))
# colnames(new.dat) <- c('corr.det', 'sel.str', 'recomload')
# new.dat <- as.data.frame(new.dat)
new.dat <- read.csv('InfinitePopRecomloadData.csv',as.is = T)
p <- ggplot(new.dat) +
  aes(x = sel.str, y = corr.det, z = recomload, fill = recomload) +
  geom_tile() + 
  # coord_equal() +
  theme(aspect.ratio=1) +
  geom_contour(color = 'white', alpha = .5) +
  #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  theme_bw() +
  labs(fill = 'Recombination \nLoad', x = 'Selection Strength on SAL', 
       y = '% Correct Sex Det.', 
       title = 'Recombination load') +
  scale_fill_gradientn(colours=viridis(256),
                       breaks=c(min(new.dat$recomload),max(new.dat$recomload)),
                       labels=c(round(min(new.dat$recomload), digits = 5),
                                round(max(new.dat$recomload), digits = 5)),
                       limits=c(min(new.dat$recomload),max(new.dat$recomload)))
p
# r<-max(abs(ggplot_build(p)$panel$ranges[[1]]$x.range))
# s<-max(abs(ggplot_build(p)$panel$ranges[[1]]$y.range))
# t<-round(max(r,s),1)
# p <- p +coord_equal(xlim=c(-t,t),ylim=c(-t,t))




# this will be the developmental load
# new.dat <-cbind(rep(seq(0,1, length.out = size), size), 
#                 rep(seq(.1,.6, length.out = size), each = size),
#                 as.vector(devmatrix))
# colnames(new.dat) <- c('dom.fac', 'sel.str', 'dev')
# new.dat <- as.data.frame(new.dat)
new.dat <- read.csv('NoisySexData.csv',as.is = T)
ggplot(new.dat) +
  aes(x = sel.str, y = dom.fac, z = dev, fill = dev) +
  geom_tile() + 
  theme(aspect.ratio=1) +
  # coord_equal() +
  geom_contour(color = 'white', alpha = .5) +
  #scale_fill_distiller(palette = 'Spectral', na.value = 'white') +
  theme_bw() +
  labs(fill = 'Noisy Sex \nLoad', x = 'Selection Strength on SAL', 
       y = 'Male Benefit Dominance Factor', 
       title = 'Noisy Sex Load') +
  scale_fill_gradientn(colours=viridis(256),
                       breaks=c(min(new.dat$dev),max(new.dat$dev)),
                       labels=c(round(min(new.dat$dev), digits = 5),
                                round(max(new.dat$dev), digits = 5)),
                       limits=c(min(new.dat$dev),max(new.dat$dev)))

