ptm <- proc.time()
N <- 10000
haps1 <- c("YGAR", "XC12", "XC13", "YC14", "YC15", "YC16")
haps2 <- c("XGAR", "XF12", "XF13", "XF14", "XF15", "XF16")
probs1 <- c(.1,.1,.1,.1,.3,.3)
probs2 <- c(.1,.3,.1,.1,.3,.1)
all.genos <- vector(length=N)
for(i in 1:10000){
  hap1 <- sample(haps1, prob=probs1, 1)
  hap2 <- sample(haps2, prob=probs2, 1)
  geno <- paste(hap1,hap2,sep=" ")
  genosex <- paste(sort(c(substr(hap1, start=1,stop=1),
                          substr(hap2, start=1,stop=1))), collapse="")
  if(genosex == "XY"){
    if(runif(1,0,1)>.9){
      geno <- paste("D", geno, collapse="")
    }else{
      geno <- paste("S", geno, collapse="")
    }
  }
  if(genosex == "XX"){
    if(runif(1,0,1)>.9){
      geno <- paste("S", geno, collapse="")
    }else{
      geno <- paste("D", geno, collapse="")
    }
  }
  all.genos[i] <- geno
}
vers.loop <- proc.time()-ptm

ptm <- proc.time()
N <- 10000
haps1 <- c("YGAR", "XC12", "XC13", "YC14", "YC15", "YC16")
haps2 <- c("XGAR", "XF12", "XF13", "XF14", "XF15", "XF16")
probs1 <- c(.1,.1,.1,.1,.3,.3)
probs2 <- c(.1,.3,.1,.1,.3,.1)
genos <- vector(length=N)all.genos <- c()
hap1 <- sample(haps1, prob=probs1, size=N, replace=T)
hap2 <- sample(haps2, prob=probs2, size=N, replace=T)
geno <- paste(hap1,hap2,sep=" ")
genosex <- paste(substr(hap1, start=1,stop=1),
                 substr(hap2, start=1,stop=1), sep=" ")
baddev <- runif(10000,0,1)>.95
phenosex <- genosex
# males = true
phenosex[genosex %in% c("X Y", "Y X", "Y Y")] <- T
phenosex[genosex == "X X"] <- F
phenosex <- as.logical(phenosex)
phenosex[baddev] <- !phenosex[baddev]
phenosex[phenosex==T] <- "S"
phenosex[phenosex=="FALSE"] <- "D"
geno <- paste(phenosex, geno)
vers.vect <- proc.time()-ptm

vers.loop
vers.vect


x <- as.data.frame(table(geno))










