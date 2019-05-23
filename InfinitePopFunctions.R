# this function runs gen for a set number of iter or until stable
# returns the iter count and the end genotype
# geno should be a vector of all of the Current allele frequencies in the order below in geno assignment
# par should be of the form c(n1,n2,h1,h2,h3,s,t,rec1,rec2,rec3,rf,gs,u)
# iter should be desired # of generations or to run until stable
# prnt is whether to print iter to console
# tol is the tolerance for what is deemed stable

NewGeno <- function(geno, par, iter = 'stable', prnt = T, tol = .000001){
  if(length(geno) != 32){stop('the genotype vector must be 32 allele freq long')}
  #if(sum(geno, na.rm = T) != 1){stop('The allele freq must sum to 1')} ****this seems to be wrong
  if(length(par) != 13){stop('The param vector must be 13 units long')}
  if(prnt != F){
    if(prnt != T){stop('Print must be true or false')}
  }
  if (!length(tol)>0 && !is.numeric(tol)){stop('tol must be a number')}
  
  # making the param vector the form to feed into gen
  param <- param.assign(par)
  
  # if iter != stable then we calculate the given number of generations
  if(iter == 'fix'){
    iter.count <- 0
    done <- F
    while(done == F){
      sec <- seq(2,32,by = 2)
      last.geno <- geno
      geno <- gen(geno,param)
      iter.count <- iter.count + 1
      if(prnt == T){cat(iter.count, '. ')}
      if(sum(geno[sec]) < tol | sum(geno[sec]) > 1-tol){
        done <- T
        fx <- T
      }else if(sum(abs(geno-last.geno)) < tol){
        done <- T
        fx <- F
      }
    }
    stop(return(list(geno,iter.count,fx)))
    # if iter = stable we run unitil the genotype reaches equilibrium
  }else if(iter == 'stable'){
    iter.count <- 0
    done <- F
    while(done == F){
      last.geno <- geno
      geno <- gen(geno,param)
      iter.count <- iter.count + 1
      if(prnt == T){cat(iter.count,' ')}
      if(sum(abs(geno-last.geno)) < tol){done <- T}
    }
  }else {
    iter.count <- iter
    for(i in 1:iter){
      if(prnt == T){cat(i,' ')}
      geno <- gen(geno,param)
    }
  }
  return(list(geno, iter.count))
}

# this function takes in the normal parameters and returns them in a form useful for the gen function
# par should be of the form c(n1,n2,h1,h2,h3,s,t,rec1,rec2,rec3,rf,gs,u)
# c(      n1,n2,h1,h2,h3,s,t,rec1,rec2,rec3,rf,gs,u)
# par <- c(.1,.1, 1, 1, 1, 0,0,0,   0,    0,  1, 0, 1)

param.assign <- function(par){
  n1 <- par[1]
  n1.c <- 1-n1
  
  n2 <- par[2]
  n2.c <- 1-n2
  
  h1 <- par[3]
  h2 <- par[4]
  h3 <- par[5]
  s <- par[6]
  t <- par[7]
  rec1 <- par[8]
  rec2 <- par[9]
  rec3 <- par[10]
  rf <- par[11]
  gs <- par[12]
  u <- par[13] # not supported
  
  w1 <- 1+s
  w2 <- 1 + h1*s
  w3 <- 1
  w4 <- 1/(1 + s)
  w5 <- 1/(1 + h1*s)
  w6 <- 1
  
  m1 <- 1
  m2 <- 1/(1 + h3*gs)
  m3 <- 1/(1 + gs)
  
  t1 <- 1
  t2 <- 1 - t
  t3 <- (1 - t)^2
  
  r1 <- rec1
  r1.c <- 1-r1
  
  r2 <- rec2
  r2.c <- 1-r2
  
  r3 <- rec3
  r3.c <- 1-r3
  
  r4 <- rec1 - rec1 * h2 * rf
  r4.c <- 1-r4
  
  r5 <- rec2 - rec2 * h2 * rf
  r5.c <- 1-r5
  
  r6 <- rec3 - rec3 * h2 * rf
  r6.c <- 1-r6
  
  r7 <- rec1 - rec1 * rf 
  r7.c <- 1-r7
  
  r8 <- rec2 - rec2 * rf
  r8.c <- 1-r8
  
  r9 <- rec3 - rec3 * rf
  r9.c <- 1-r9
  
  param <- c(n1,n1.c,n2,n2.c,w1,w2,w3,w4,w5,w6,m1,m2,m3,t1,t2,t3,r1,r1.c,r2,r2.c,r3,r3.c,r4,r4.c,
             r5,r5.c,r6,r6.c,r7,r7.c,r8,r8.c,r9,r9.c)
  return(param)
}

# this function takes a vector of allele freq and parameters from param.assign and returns the 
# allele frequencies of the following generation
# geno should be a vector of all of the Current allele frequencies in the order below in geno assignment
# param should be of the form  c(n1,n1.c,n2,n2.c,w1,w2,w3,w4,w5,w6,m1,m2,m3,t1,t2,t3,r1,r1.c,r2,r2.c,
# r3,r3.c,r4,r4.c,r5,r5.c,r6,r6.c,r7,r7.c,r8,r8.c,r9,r9.c)
# this is how it is returned from param.assign
# if return.fit = T this function returns a vector of the male and female mean fitness

gen <- function(geno,param, return.fit = F, return.abs.fit = F){
  ##### geno assignment #####
  XGARs <- geno[1]
  XGArs <- geno[2]
  XGaRs <- geno[3]
  XGars <- geno[4]
  XgARs <- geno[5]
  XgArs <- geno[6]
  XgaRs <- geno[7]
  Xgars <- geno[8]
  YGARs <- geno[9]
  YGArs <- geno[10]
  YGaRs <- geno[11]
  YGars <- geno[12]
  YgARs <- geno[13]
  YgArs <- geno[14]
  YgaRs <- geno[15]
  Ygars <- geno[16]
  XGARe <- geno[17]
  XGAre <- geno[18]
  XGaRe <- geno[19]
  XGare <- geno[20]
  XgARe <- geno[21]
  XgAre <- geno[22]
  XgaRe <- geno[23]
  Xgare <- geno[24]
  YGARe <- geno[25]
  YGAre <- geno[26]
  YGaRe <- geno[27]
  YGare <- geno[28]
  YgARe <- geno[29]
  YgAre <- geno[30]
  YgaRe <- geno[31]
  Ygare <- geno[32]
  
  ##### param assignment #####
  n1 <- param[1]
  n1.c <- param[2]
  n2 <- param[3]
  n2.c <- param[4]
  w1 <- param[5]
  w2 <- param[6]
  w3 <- param[7]
  w4 <- param[8]
  w5 <- param[9]
  w6 <- param[10]
  m1 <- param[11]
  m2 <- param[12]
  m3 <- param[13]
  t1 <- param[14]
  t2 <- param[15]
  t3 <- param[16]
  r1 <- param[17]
  r1.c <- param[18]
  r2 <- param[19]
  r2.c <- param[20]
  r3 <- param[21]
  r3.c <- param[22]
  r4 <- param[23]
  r4.c <- param[24]
  r5 <- param[25]
  r5.c <- param[26]
  r6 <- param[27]
  r6.c <- param[28]
  r7 <- param[29]
  r7.c <- param[30]
  r8 <- param[31]
  r8.c <- param[32]
  r9 <- param[33]
  r9.c <- param[34]
  
  ##### male.fit.bar #####
  male.fit.bar <- 2 * (XGARs*XGARe * n2.c * w1 * m1 * t2 +
                         XGARs*XGAre * n2.c * w1 * m1 * t2 +
                         XGARs*XGaRe * n2.c * w2 * m1 * t2 +
                         XGARs*XGare * n2.c * w2 * m1 * t2 +
                         XGARs*XgARe * n2.c * w1 * m2 * t2 +
                         XGARs*XgAre * n2.c * w1 * m2 * t2 +
                         XGARs*XgaRe * n2.c * w2 * m2 * t2 +
                         XGARs*Xgare * n2.c * w2 * m2 * t2 +
                         XGARs*YGARe * n1 * w1 * m1 * t1 +
                         XGARs*YGAre * n1 * w1 * m1 * t1 +
                         XGARs*YGaRe * n1 * w2 * m1 * t1 +
                         XGARs*YGare * n1 * w2 * m1 * t1 +
                         XGARs*YgARe * n1 * w1 * m2 * t1 +
                         XGARs*YgAre * n1 * w1 * m2 * t1 +
                         XGARs*YgaRe * n1 * w2 * m2 * t1 +
                         XGARs*Ygare * n1 * w2 * m2 * t1 +
                         # row 2
                         XGArs*XGARe * n2.c * w1 * m1 * t2 +
                         XGArs*XGAre * n2.c * w1 * m1 * t2 +
                         XGArs*XGaRe * n2.c * w2 * m1 * t2 +
                         XGArs*XGare * n2.c * w2 * m1 * t2 +
                         XGArs*XgARe * n2.c * w1 * m2 * t2 +
                         XGArs*XgAre * n2.c * w1 * m2 * t2 +
                         XGArs*XgaRe * n2.c * w2 * m2 * t2 +
                         XGArs*Xgare * n2.c * w2 * m2 * t2 +
                         XGArs*YGARe * n1 * w1 * m1 * t1 +
                         XGArs*YGAre * n1 * w1 * m1 * t1 +
                         XGArs*YGaRe * n1 * w2 * m1 * t1 +
                         XGArs*YGare * n1 * w2 * m1 * t1 +
                         XGArs*YgARe * n1 * w1 * m2 * t1 +
                         XGArs*YgAre * n1 * w1 * m2 * t1 +
                         XGArs*YgaRe * n1 * w2 * m2 * t1 +
                         XGArs*Ygare * n1 * w2 * m2 * t1 +
                         # row 3
                         XGaRs*XGARe * n2.c * w2 * m1 * t2 +
                         XGaRs*XGAre * n2.c * w2 * m1 * t2 +
                         XGaRs*XGaRe * n2.c * w3 * m1 * t2 +
                         XGaRs*XGare * n2.c * w3 * m1 * t2 +
                         XGaRs*XgARe * n2.c * w2 * m2 * t2 +
                         XGaRs*XgAre * n2.c * w2 * m2 * t2 +
                         XGaRs*XgaRe * n2.c * w3 * m2 * t2 +
                         XGaRs*Xgare * n2.c * w3 * m2 * t2 +
                         XGaRs*YGARe * n1 * w2 * m1 * t1 +
                         XGaRs*YGAre * n1 * w2 * m1 * t1 +
                         XGaRs*YGaRe * n1 * w3 * m1 * t1 +
                         XGaRs*YGare * n1 * w3 * m1 * t1 +
                         XGaRs*YgARe * n1 * w2 * m2 * t1 +
                         XGaRs*YgAre * n1 * w2 * m2 * t1 +
                         XGaRs*YgaRe * n1 * w3 * m2 * t1 +
                         XGaRs*Ygare * n1 * w3 * m2 * t1 +
                         # row 4
                         XGars*XGARe * n2.c * w2 * m1 * t2 +
                         XGars*XGAre * n2.c * w2 * m1 * t2 +
                         XGars*XGaRe * n2.c * w3 * m1 * t2 +
                         XGars*XGare * n2.c * w3 * m1 * t2 +
                         XGars*XgARe * n2.c * w2 * m2 * t2 +
                         XGars*XgAre * n2.c * w2 * m2 * t2 +
                         XGars*XgaRe * n2.c * w3 * m2 * t2 +
                         XGars*Xgare * n2.c * w3 * m2 * t2 +
                         XGars*YGARe * n1 * w2 * m1 * t1 +
                         XGars*YGAre * n1 * w2 * m1 * t1 +
                         XGars*YGaRe * n1 * w3 * m1 * t1 +
                         XGars*YGare * n1 * w3 * m1 * t1 +
                         XGars*YgARe * n1 * w2 * m2 * t1 +
                         XGars*YgAre * n1 * w2 * m2 * t1 +
                         XGars*YgaRe * n1 * w3 * m2 * t1 +
                         XGars*Ygare * n1 * w3 * m2 * t1 +
                         # row 5
                         XgARs*XGARe * n2.c * w1 * m2 * t2 +
                         XgARs*XGAre * n2.c * w1 * m2 * t2 +
                         XgARs*XGaRe * n2.c * w2 * m2 * t2 +
                         XgARs*XGare * n2.c * w2 * m2 * t2 +
                         XgARs*XgARe * n2.c * w1 * m3 * t2 +
                         XgARs*XgAre * n2.c * w1 * m3 * t2 +
                         XgARs*XgaRe * n2.c * w2 * m3 * t2 +
                         XgARs*Xgare * n2.c * w2 * m3 * t2 +
                         XgARs*YGARe * n1 * w1 * m2 * t1 +
                         XgARs*YGAre * n1 * w1 * m2 * t1 +
                         XgARs*YGaRe * n1 * w2 * m2 * t1 +
                         XgARs*YGare * n1 * w2 * m2 * t1 +
                         XgARs*YgARe * n1 * w1 * m3 * t1 +
                         XgARs*YgAre * n1 * w1 * m3 * t1 +
                         XgARs*YgaRe * n1 * w2 * m3 * t1 +
                         XgARs*Ygare * n1 * w2 * m3 * t1 +
                         # row 6
                         XgArs*XGARe * n2.c * w1 * m2 * t2 +
                         XgArs*XGAre * n2.c * w1 * m2 * t2 +
                         XgArs*XGaRe * n2.c * w2 * m2 * t2 +
                         XgArs*XGare * n2.c * w2 * m2 * t2 +
                         XgArs*XgARe * n2.c * w1 * m3 * t2 +
                         XgArs*XgAre * n2.c * w1 * m3 * t2 +
                         XgArs*XgaRe * n2.c * w2 * m3 * t2 +
                         XgArs*Xgare * n2.c * w2 * m3 * t2 +
                         XgArs*YGARe * n1 * w1 * m2 * t1 +
                         XgArs*YGAre * n1 * w1 * m2 * t1 +
                         XgArs*YGaRe * n1 * w2 * m2 * t1 +
                         XgArs*YGare * n1 * w2 * m2 * t1 +
                         XgArs*YgARe * n1 * w1 * m3 * t1 +
                         XgArs*YgAre * n1 * w1 * m3 * t1 +
                         XgArs*YgaRe * n1 * w2 * m3 * t1 +
                         XgArs*Ygare * n1 * w2 * m3 * t1 +
                         # row 7
                         XgaRs*XGARe * n2.c * w2 * m2 * t2 +
                         XgaRs*XGAre * n2.c * w2 * m2 * t2 +
                         XgaRs*XGaRe * n2.c * w3 * m2 * t2 +
                         XgaRs*XGare * n2.c * w3 * m2 * t2 +
                         XgaRs*XgARe * n2.c * w2 * m3 * t2 +
                         XgaRs*XgAre * n2.c * w2 * m3 * t2 +
                         XgaRs*XgaRe * n2.c * w3 * m3 * t2 +
                         XgaRs*Xgare * n2.c * w3 * m3 * t2 +
                         XgaRs*YGARe * n1 * w2 * m2 * t1 +
                         XgaRs*YGAre * n1 * w2 * m2 * t1 +
                         XgaRs*YGaRe * n1 * w3 * m2 * t1 +
                         XgaRs*YGare * n1 * w3 * m2 * t1 +
                         XgaRs*YgARe * n1 * w2 * m3 * t1 +
                         XgaRs*YgAre * n1 * w2 * m3 * t1 +
                         XgaRs*YgaRe * n1 * w3 * m3 * t1 +
                         XgaRs*Ygare * n1 * w3 * m3 * t1 +
                         # row 8
                         Xgars*XGARe * n2.c * w2 * m2 * t2 +
                         Xgars*XGAre * n2.c * w2 * m2 * t2 +
                         Xgars*XGaRe * n2.c * w3 * m2 * t2 +
                         Xgars*XGare * n2.c * w3 * m2 * t2 +
                         Xgars*XgARe * n2.c * w2 * m3 * t2 +
                         Xgars*XgAre * n2.c * w2 * m3 * t2 +
                         Xgars*XgaRe * n2.c * w3 * m3 * t2 +
                         Xgars*Xgare * n2.c * w3 * m3 * t2 +
                         Xgars*YGARe * n1 * w2 * m2 * t1 +
                         Xgars*YGAre * n1 * w2 * m2 * t1 +
                         Xgars*YGaRe * n1 * w3 * m2 * t1 +
                         Xgars*YGare * n1 * w3 * m2 * t1 +
                         Xgars*YgARe * n1 * w2 * m3 * t1 +
                         Xgars*YgAre * n1 * w2 * m3 * t1 +
                         Xgars*YgaRe * n1 * w3 * m3 * t1 +
                         Xgars*Ygare * n1 * w3 * m3 * t1 +
                         # row 9
                         YGARs*XGARe * n1 * w1 * m1 * t1 +
                         YGARs*XGAre * n1 * w1 * m1 * t1 +
                         YGARs*XGaRe * n1 * w2 * m1 * t1 +
                         YGARs*XGare * n1 * w2 * m1 * t1 +
                         YGARs*XgARe * n1 * w1 * m2 * t1 +
                         YGARs*XgAre * n1 * w1 * m2 * t1 +
                         YGARs*XgaRe * n1 * w2 * m2 * t1 +
                         YGARs*Xgare * n1 * w2 * m2 * t1 +
                         YGARs*YGARe * n1 * w1 * m1 * t2 +
                         YGARs*YGAre * n1 * w1 * m1 * t2 +
                         YGARs*YGaRe * n1 * w2 * m1 * t2 +
                         YGARs*YGare * n1 * w2 * m1 * t2 +
                         YGARs*YgARe * n1 * w1 * m2 * t2 +
                         YGARs*YgAre * n1 * w1 * m2 * t2 +
                         YGARs*YgaRe * n1 * w2 * m2 * t2 +
                         YGARs*Ygare * n1 * w2 * m2 * t2 +
                         # row 10
                         YGArs*XGARe * n1 * w1 * m1 * t1 +
                         YGArs*XGAre * n1 * w1 * m1 * t1 +
                         YGArs*XGaRe * n1 * w2 * m1 * t1 +
                         YGArs*XGare * n1 * w2 * m1 * t1 +
                         YGArs*XgARe * n1 * w1 * m2 * t1 +
                         YGArs*XgAre * n1 * w1 * m2 * t1 +
                         YGArs*XgaRe * n1 * w2 * m2 * t1 +
                         YGArs*Xgare * n1 * w2 * m2 * t1 +
                         YGArs*YGARe * n1 * w1 * m1 * t2 +
                         YGArs*YGAre * n1 * w1 * m1 * t2 +
                         YGArs*YGaRe * n1 * w2 * m1 * t2 +
                         YGArs*YGare * n1 * w2 * m1 * t2 +
                         YGArs*YgARe * n1 * w1 * m2 * t2 +
                         YGArs*YgAre * n1 * w1 * m2 * t2 +
                         YGArs*YgaRe * n1 * w2 * m2 * t2 +
                         YGArs*Ygare * n1 * w2 * m2 * t2 +
                         # row 11
                         YGaRs*XGARe * n1 * w2 * m1 * t1 +
                         YGaRs*XGAre * n1 * w2 * m1 * t1 +
                         YGaRs*XGaRe * n1 * w3 * m1 * t1 +
                         YGaRs*XGare * n1 * w3 * m1 * t1 +
                         YGaRs*XgARe * n1 * w2 * m2 * t1 +
                         YGaRs*XgAre * n1 * w2 * m2 * t1 +
                         YGaRs*XgaRe * n1 * w3 * m2 * t1 +
                         YGaRs*Xgare * n1 * w3 * m2 * t1 +
                         YGaRs*YGARe * n1 * w2 * m1 * t2 +
                         YGaRs*YGAre * n1 * w2 * m1 * t2 +
                         YGaRs*YGaRe * n1 * w3 * m1 * t2 +
                         YGaRs*YGare * n1 * w3 * m1 * t2 +
                         YGaRs*YgARe * n1 * w2 * m2 * t2 +
                         YGaRs*YgAre * n1 * w2 * m2 * t2 +
                         YGaRs*YgaRe * n1 * w3 * m2 * t2 +
                         YGaRs*Ygare * n1 * w3 * m2 * t2 +
                         # row 12
                         YGars*XGARe * n1 * w2 * m1 * t1 +
                         YGars*XGAre * n1 * w2 * m1 * t1 +
                         YGars*XGaRe * n1 * w3 * m1 * t1 +
                         YGars*XGare * n1 * w3 * m1 * t1 +
                         YGars*XgARe * n1 * w2 * m2 * t1 +
                         YGars*XgAre * n1 * w2 * m2 * t1 +
                         YGars*XgaRe * n1 * w3 * m2 * t1 +
                         YGars*Xgare * n1 * w3 * m2 * t1 +
                         YGars*YGARe * n1 * w2 * m1 * t2 +
                         YGars*YGAre * n1 * w2 * m1 * t2 +
                         YGars*YGaRe * n1 * w3 * m1 * t2 +
                         YGars*YGare * n1 * w3 * m1 * t2 +
                         YGars*YgARe * n1 * w2 * m2 * t2 +
                         YGars*YgAre * n1 * w2 * m2 * t2 +
                         YGars*YgaRe * n1 * w3 * m2 * t2 +
                         YGars*Ygare * n1 * w3 * m2 * t2 +
                         # row 13
                         YgARs*XGARe * n1 * w1 * m2 * t1 +
                         YgARs*XGAre * n1 * w1 * m2 * t1 +
                         YgARs*XGaRe * n1 * w2 * m2 * t1 +
                         YgARs*XGare * n1 * w2 * m2 * t1 +
                         YgARs*XgARe * n1 * w1 * m3 * t1 +
                         YgARs*XgAre * n1 * w1 * m3 * t1 +
                         YgARs*XgaRe * n1 * w2 * m3 * t1 +
                         YgARs*Xgare * n1 * w2 * m3 * t1 +
                         YgARs*YGARe * n1 * w1 * m2 * t2 +
                         YgARs*YGAre * n1 * w1 * m2 * t2 +
                         YgARs*YGaRe * n1 * w2 * m2 * t2 +
                         YgARs*YGare * n1 * w2 * m2 * t2 +
                         YgARs*YgARe * n1 * w1 * m3 * t2 +
                         YgARs*YgAre * n1 * w1 * m3 * t2 +
                         YgARs*YgaRe * n1 * w2 * m3 * t2 +
                         YgARs*Ygare * n1 * w2 * m3 * t2 +
                         # row 14
                         YgArs*XGARe * n1 * w1 * m2 * t1 +
                         YgArs*XGAre * n1 * w1 * m2 * t1 +
                         YgArs*XGaRe * n1 * w2 * m2 * t1 +
                         YgArs*XGare * n1 * w2 * m2 * t1 +
                         YgArs*XgARe * n1 * w1 * m3 * t1 +
                         YgArs*XgAre * n1 * w1 * m3 * t1 +
                         YgArs*XgaRe * n1 * w2 * m3 * t1 +
                         YgArs*Xgare * n1 * w2 * m3 * t1 +
                         YgArs*YGARe * n1 * w1 * m2 * t2 +
                         YgArs*YGAre * n1 * w1 * m2 * t2 +
                         YgArs*YGaRe * n1 * w2 * m2 * t2 +
                         YgArs*YGare * n1 * w2 * m2 * t2 +
                         YgArs*YgARe * n1 * w1 * m3 * t2 +
                         YgArs*YgAre * n1 * w1 * m3 * t2 +
                         YgArs*YgaRe * n1 * w2 * m3 * t2 +
                         YgArs*Ygare * n1 * w2 * m3 * t2 +
                         # row 15
                         YgaRs*XGARe * n1 * w2 * m2 * t1 +
                         YgaRs*XGAre * n1 * w2 * m2 * t1 +
                         YgaRs*XGaRe * n1 * w3 * m2 * t1 +
                         YgaRs*XGare * n1 * w3 * m2 * t1 +
                         YgaRs*XgARe * n1 * w2 * m3 * t1 +
                         YgaRs*XgAre * n1 * w2 * m3 * t1 +
                         YgaRs*XgaRe * n1 * w3 * m3 * t1 +
                         YgaRs*Xgare * n1 * w3 * m3 * t1 +
                         YgaRs*YGARe * n1 * w2 * m2 * t2 +
                         YgaRs*YGAre * n1 * w2 * m2 * t2 +
                         YgaRs*YGaRe * n1 * w3 * m2 * t2 +
                         YgaRs*YGare * n1 * w3 * m2 * t2 +
                         YgaRs*YgARe * n1 * w2 * m3 * t2 +
                         YgaRs*YgAre * n1 * w2 * m3 * t2 +
                         YgaRs*YgaRe * n1 * w3 * m3 * t2 +
                         YgaRs*Ygare * n1 * w3 * m3 * t2 +
                         # row 16
                         Ygars*XGARe * n1 * w2 * m2 * t1 +
                         Ygars*XGAre * n1 * w2 * m2 * t1 +
                         Ygars*XGaRe * n1 * w3 * m2 * t1 +
                         Ygars*XGare * n1 * w3 * m2 * t1 +
                         Ygars*XgARe * n1 * w2 * m3 * t1 +
                         Ygars*XgAre * n1 * w2 * m3 * t1 +
                         Ygars*XgaRe * n1 * w3 * m3 * t1 +
                         Ygars*Xgare * n1 * w3 * m3 * t1 +
                         Ygars*YGARe * n1 * w2 * m2 * t2 +
                         Ygars*YGAre * n1 * w2 * m2 * t2 +
                         Ygars*YGaRe * n1 * w3 * m2 * t2 +
                         Ygars*YGare * n1 * w3 * m2 * t2 +
                         Ygars*YgARe * n1 * w2 * m3 * t2 +
                         Ygars*YgAre * n1 * w2 * m3 * t2 +
                         Ygars*YgaRe * n1 * w3 * m3 * t2 +
                         Ygars*Ygare * n1 * w3 * m3 * t2)
  
  ##### female.fit.bar #####
  female.fit.bar <- 2 * (XGARs*XGARe * n2 * w4 * m1 * t1 +
                           XGARs*XGAre * n2 * w4 * m1 * t1 +
                           XGARs*XGaRe * n2 * w5 * m1 * t1 +
                           XGARs*XGare * n2 * w5 * m1 * t1 +
                           XGARs*XgARe * n2 * w4 * m2 * t1 +
                           XGARs*XgAre * n2 * w4 * m2 * t1 +
                           XGARs*XgaRe * n2 * w5 * m2 * t1 +
                           XGARs*Xgare * n2 * w5 * m2 * t1 +
                           XGARs*YGARe * n1.c * w4 * m1 * t2 +
                           XGARs*YGAre * n1.c * w4 * m1 * t2 +
                           XGARs*YGaRe * n1.c * w5 * m1 * t2 +
                           XGARs*YGare * n1.c * w5 * m1 * t2 +
                           XGARs*YgARe * n1.c * w4 * m2 * t2 +
                           XGARs*YgAre * n1.c * w4 * m2 * t2 +
                           XGARs*YgaRe * n1.c * w5 * m2 * t2 +
                           XGARs*Ygare * n1.c * w5 * m2 * t2 +
                           # row 2
                           XGArs*XGARe * n2 * w4 * m1 * t1 +
                           XGArs*XGAre * n2 * w4 * m1 * t1 +
                           XGArs*XGaRe * n2 * w5 * m1 * t1 +
                           XGArs*XGare * n2 * w5 * m1 * t1 +
                           XGArs*XgARe * n2 * w4 * m2 * t1 +
                           XGArs*XgAre * n2 * w4 * m2 * t1 +
                           XGArs*XgaRe * n2 * w5 * m2 * t1 +
                           XGArs*Xgare * n2 * w5 * m2 * t1 +
                           XGArs*YGARe * n1.c * w4 * m1 * t2 +
                           XGArs*YGAre * n1.c * w4 * m1 * t2 +
                           XGArs*YGaRe * n1.c * w5 * m1 * t2 +
                           XGArs*YGare * n1.c * w5 * m1 * t2 +
                           XGArs*YgARe * n1.c * w4 * m2 * t2 +
                           XGArs*YgAre * n1.c * w4 * m2 * t2 +
                           XGArs*YgaRe * n1.c * w5 * m2 * t2 +
                           XGArs*Ygare * n1.c * w5 * m2 * t2 +
                           # row 3
                           XGaRs*XGARe * n2 * w5 * m1 * t1 +
                           XGaRs*XGAre * n2 * w5 * m1 * t1 +
                           XGaRs*XGaRe * n2 * w6 * m1 * t1 +
                           XGaRs*XGare * n2 * w6 * m1 * t1 +
                           XGaRs*XgARe * n2 * w5 * m2 * t1 +
                           XGaRs*XgAre * n2 * w5 * m2 * t1 +
                           XGaRs*XgaRe * n2 * w6 * m2 * t1 +
                           XGaRs*Xgare * n2 * w6 * m2 * t1 +
                           XGaRs*YGARe * n1.c * w5 * m1 * t2 +
                           XGaRs*YGAre * n1.c * w5 * m1 * t2 +
                           XGaRs*YGaRe * n1.c * w6 * m1 * t2 +
                           XGaRs*YGare * n1.c * w6 * m1 * t2 +
                           XGaRs*YgARe * n1.c * w5 * m2 * t2 +
                           XGaRs*YgAre * n1.c * w5 * m2 * t2 +
                           XGaRs*YgaRe * n1.c * w6 * m2 * t2 + 
                           XGaRs*Ygare * n1.c * w6 * m2 * t2 +
                           # row 4
                           XGars*XGARe * n2 * w5 * m1 * t1 +
                           XGars*XGAre * n2 * w5 * m1 * t1 +
                           XGars*XGaRe * n2 * w6 * m1 * t1 +
                           XGars*XGare * n2 * w6 * m1 * t1 +
                           XGars*XgARe * n2 * w5 * m2 * t1 +
                           XGars*XgAre * n2 * w5 * m2 * t1 +
                           XGars*XgaRe * n2 * w6 * m2 * t1 +
                           XGars*Xgare * n2 * w6 * m2 * t1 +
                           XGars*YGARe * n1.c * w5 * m1 * t2 +
                           XGars*YGAre * n1.c * w5 * m1 * t2 +
                           XGars*YGaRe * n1.c * w6 * m1 * t2 +
                           XGars*YGare * n1.c * w6 * m1 * t2 +
                           XGars*YgARe * n1.c * w5 * m2 * t2 +
                           XGars*YgAre * n1.c * w5 * m2 * t2 +
                           XGars*YgaRe * n1.c * w6 * m2 * t2 +
                           XGars*Ygare * n1.c * w6 * m2 * t2 +
                           # row 5
                           XgARs*XGARe * n2 * w4 * m2 * t1 +
                           XgARs*XGAre * n2 * w4 * m2 * t1 +
                           XgARs*XGaRe * n2 * w5 * m2 * t1 +
                           XgARs*XGare * n2 * w5 * m2 * t1 +
                           XgARs*XgARe * n2 * w4 * m3 * t1 +
                           XgARs*XgAre * n2 * w4 * m3 * t1 +
                           XgARs*XgaRe * n2 * w5 * m3 * t1 +
                           XgARs*Xgare * n2 * w5 * m3 * t1 +
                           XgARs*YGARe * n1.c * w4 * m2 * t2 +
                           XgARs*YGAre * n1.c * w4 * m2 * t2 +
                           XgARs*YGaRe * n1.c * w5 * m2 * t2 +
                           XgARs*YGare * n1.c * w5 * m2 * t2 +
                           XgARs*YgARe * n1.c * w4 * m3 * t2 +
                           XgARs*YgAre * n1.c * w4 * m3 * t2 +
                           XgARs*YgaRe * n1.c * w5 * m3 * t2 +
                           XgARs*Ygare * n1.c * w5 * m3 * t2 +
                           # row 6
                           XgArs*XGARe * n2 * w4 * m2 * t1 +
                           XgArs*XGAre * n2 * w4 * m2 * t1 +
                           XgArs*XGaRe * n2 * w5 * m2 * t1 +
                           XgArs*XGare * n2 * w5 * m2 * t1 +
                           XgArs*XgARe * n2 * w4 * m3 * t1 +
                           XgArs*XgAre * n2 * w4 * m3 * t1 +
                           XgArs*XgaRe * n2 * w5 * m3 * t1 +
                           XgArs*Xgare * n2 * w5 * m3 * t1 +
                           XgArs*YGARe * n1.c * w4 * m2 * t2 +
                           XgArs*YGAre * n1.c * w4 * m2 * t2 +
                           XgArs*YGaRe * n1.c * w5 * m2 * t2 +
                           XgArs*YGare * n1.c * w5 * m2 * t2 +
                           XgArs*YgARe * n1.c * w4 * m3 * t2 +
                           XgArs*YgAre * n1.c * w4 * m3 * t2 +
                           XgArs*YgaRe * n1.c * w5 * m3 * t2 +
                           XgArs*Ygare * n1.c * w5 * m3 * t2 +
                           # row 7
                           XgaRs*XGARe * n2 * w5 * m2 * t1 +
                           XgaRs*XGAre * n2 * w5 * m2 * t1 +
                           XgaRs*XGaRe * n2 * w6 * m2 * t1 +
                           XgaRs*XGare * n2 * w6 * m2 * t1 +
                           XgaRs*XgARe * n2 * w5 * m3 * t1 +
                           XgaRs*XgAre * n2 * w5 * m3 * t1 +
                           XgaRs*XgaRe * n2 * w6 * m3 * t1 +
                           XgaRs*Xgare * n2 * w6 * m3 * t1 +
                           XgaRs*YGARe * n1.c * w5 * m2 * t2 +
                           XgaRs*YGAre * n1.c * w5 * m2 * t2 +
                           XgaRs*YGaRe * n1.c * w6 * m2 * t2 +
                           XgaRs*YGare * n1.c * w6 * m2 * t2 +
                           XgaRs*YgARe * n1.c * w5 * m3 * t2 +
                           XgaRs*YgAre * n1.c * w5 * m3 * t2 +
                           XgaRs*YgaRe * n1.c * w6 * m3 * t2 +
                           XgaRs*Ygare * n1.c * w6 * m3 * t2 +
                           # row 8
                           Xgars*XGARe * n2 * w5 * m2 * t1 +
                           Xgars*XGAre * n2 * w5 * m2 * t1 +
                           Xgars*XGaRe * n2 * w6 * m2 * t1 +
                           Xgars*XGare * n2 * w6 * m2 * t1 +
                           Xgars*XgARe * n2 * w5 * m3 * t1 +
                           Xgars*XgAre * n2 * w5 * m3 * t1 +
                           Xgars*XgaRe * n2 * w6 * m3 * t1 +
                           Xgars*Xgare * n2 * w6 * m3 * t1 +
                           Xgars*YGARe * n1.c * w5 * m2 * t2 +
                           Xgars*YGAre * n1.c * w5 * m2 * t2 +
                           Xgars*YGaRe * n1.c * w6 * m2 * t2 +
                           Xgars*YGare * n1.c * w6 * m2 * t2 +
                           Xgars*YgARe * n1.c * w5 * m3 * t2 +
                           Xgars*YgAre * n1.c * w5 * m3 * t2 +
                           Xgars*YgaRe * n1.c * w6 * m3 * t2 +
                           Xgars*Ygare * n1.c * w6 * m3 * t2 +
                           # row 9
                           YGARs*XGARe * n1.c * w4 * m1 * t2 +
                           YGARs*XGAre * n1.c * w4 * m1 * t2 +
                           YGARs*XGaRe * n1.c * w5 * m1 * t2 +
                           YGARs*XGare * n1.c * w5 * m1 * t2 +
                           YGARs*XgARe * n1.c * w4 * m2 * t2 +
                           YGARs*XgAre * n1.c * w4 * m2 * t2 +
                           YGARs*XgaRe * n1.c * w5 * m2 * t2 +
                           YGARs*Xgare * n1.c * w5 * m2 * t2 +
                           YGARs*YGARe * n1.c * w4 * m1 * t3 +
                           YGARs*YGAre * n1.c * w4 * m1 * t3 +
                           YGARs*YGaRe * n1.c * w5 * m1 * t3 +
                           YGARs*YGare * n1.c * w5 * m1 * t3 +
                           YGARs*YgARe * n1.c * w4 * m2 * t3 +
                           YGARs*YgAre * n1.c * w4 * m2 * t3 +
                           YGARs*YgaRe * n1.c * w5 * m2 * t3 +
                           YGARs*Ygare * n1.c * w5 * m2 * t3 +
                           # row 10
                           YGArs*XGARe * n1.c * w4 * m1 * t2 +
                           YGArs*XGAre * n1.c * w4 * m1 * t2 +
                           YGArs*XGaRe * n1.c * w5 * m1 * t2 +
                           YGArs*XGare * n1.c * w5 * m1 * t2 +
                           YGArs*XgARe * n1.c * w4 * m2 * t2 +
                           YGArs*XgAre * n1.c * w4 * m2 * t2 +
                           YGArs*XgaRe * n1.c * w5 * m2 * t2 +
                           YGArs*Xgare * n1.c * w5 * m2 * t2 +
                           YGArs*YGARe * n1.c * w4 * m1 * t3 +
                           YGArs*YGAre * n1.c * w4 * m1 * t3 +
                           YGArs*YGaRe * n1.c * w5 * m1 * t3 +
                           YGArs*YGare * n1.c * w5 * m1 * t3 +
                           YGArs*YgARe * n1.c * w4 * m2 * t3 +
                           YGArs*YgAre * n1.c * w4 * m2 * t3 +
                           YGArs*YgaRe * n1.c * w5 * m2 * t3 +
                           YGArs*Ygare * n1.c * w5 * m2 * t3 +
                           # row 11
                           YGaRs*XGARe * n1.c * w5 * m1 * t2 +
                           YGaRs*XGAre * n1.c * w5 * m1 * t2 +
                           YGaRs*XGaRe * n1.c * w6 * m1 * t2 +
                           YGaRs*XGare * n1.c * w6 * m1 * t2 +
                           YGaRs*XgARe * n1.c * w5 * m2 * t2 +
                           YGaRs*XgAre * n1.c * w5 * m2 * t2 +
                           YGaRs*XgaRe * n1.c * w6 * m2 * t2 +
                           YGaRs*Xgare * n1.c * w6 * m2 * t2 +
                           YGaRs*YGARe * n1.c * w5 * m1 * t3 +
                           YGaRs*YGAre * n1.c * w5 * m1 * t3 +
                           YGaRs*YGaRe * n1.c * w6 * m1 * t3 +
                           YGaRs*YGare * n1.c * w6 * m1 * t3 +
                           YGaRs*YgARe * n1.c * w5 * m2 * t3 +
                           YGaRs*YgAre * n1.c * w5 * m2 * t3 +
                           YGaRs*YgaRe * n1.c * w6 * m2 * t3 +
                           YGaRs*Ygare * n1.c * w6 * m2 * t3 +
                           # row 12
                           YGars*XGARe * n1.c * w5 * m1 * t2 +
                           YGars*XGAre * n1.c * w5 * m1 * t2 +
                           YGars*XGaRe * n1.c * w6 * m1 * t2 +
                           YGars*XGare * n1.c * w6 * m1 * t2 +
                           YGars*XgARe * n1.c * w5 * m2 * t2 +
                           YGars*XgAre * n1.c * w5 * m2 * t2 +
                           YGars*XgaRe * n1.c * w6 * m2 * t2 +
                           YGars*Xgare * n1.c * w6 * m2 * t2 +
                           YGars*YGARe * n1.c * w5 * m1 * t3 +
                           YGars*YGAre * n1.c * w5 * m1 * t3 +
                           YGars*YGaRe * n1.c * w6 * m1 * t3 +
                           YGars*YGare * n1.c * w6 * m1 * t3 +
                           YGars*YgARe * n1.c * w5 * m2 * t3 +
                           YGars*YgAre * n1.c * w5 * m2 * t3 +
                           YGars*YgaRe * n1.c * w6 * m2 * t3 +
                           YGars*Ygare * n1.c * w6 * m2 * t3 +
                           # row 13
                           YgARs*XGARe * n1.c * w4 * m2 * t2 +
                           YgARs*XGAre * n1.c * w4 * m2 * t2 +
                           YgARs*XGaRe * n1.c * w5 * m2 * t2 +
                           YgARs*XGare * n1.c * w5 * m2 * t2 +
                           YgARs*XgARe * n1.c * w4 * m3 * t2 +
                           YgARs*XgAre * n1.c * w4 * m3 * t2 +
                           YgARs*XgaRe * n1.c * w5 * m3 * t2 +
                           YgARs*Xgare * n1.c * w5 * m3 * t2 +
                           YgARs*YGARe * n1.c * w4 * m2 * t3 +
                           YgARs*YGAre * n1.c * w4 * m2 * t3 +
                           YgARs*YGaRe * n1.c * w5 * m2 * t3 +
                           YgARs*YGare * n1.c * w5 * m2 * t3 +
                           YgARs*YgARe * n1.c * w4 * m3 * t3 +
                           YgARs*YgAre * n1.c * w4 * m3 * t3 +
                           YgARs*YgaRe * n1.c * w5 * m3 * t3 +
                           YgARs*Ygare * n1.c * w5 * m3 * t3 +
                           # row 14
                           YgArs*XGARe * n1.c * w4 * m2 * t2 +
                           YgArs*XGAre * n1.c * w4 * m2 * t2 +
                           YgArs*XGaRe * n1.c * w5 * m2 * t2 +
                           YgArs*XGare * n1.c * w5 * m2 * t2 +
                           YgArs*XgARe * n1.c * w4 * m3 * t2 +
                           YgArs*XgAre * n1.c * w4 * m3 * t2 +
                           YgArs*XgaRe * n1.c * w5 * m3 * t2 +
                           YgArs*Xgare * n1.c * w5 * m3 * t2 +
                           YgArs*YGARe * n1.c * w4 * m2 * t3 +
                           YgArs*YGAre * n1.c * w4 * m2 * t3 +
                           YgArs*YGaRe * n1.c * w5 * m2 * t3 +
                           YgArs*YGare * n1.c * w5 * m2 * t3 +
                           YgArs*YgARe * n1.c * w4 * m3 * t3 +
                           YgArs*YgAre * n1.c * w4 * m3 * t3 +
                           YgArs*YgaRe * n1.c * w5 * m3 * t3 +
                           YgArs*Ygare * n1.c * w5 * m3 * t3 +
                           # row 15
                           YgaRs*XGARe * n1.c * w5 * m2 * t2 +
                           YgaRs*XGAre * n1.c * w5 * m2 * t2 +
                           YgaRs*XGaRe * n1.c * w6 * m2 * t2 +
                           YgaRs*XGare * n1.c * w6 * m2 * t2 +
                           YgaRs*XgARe * n1.c * w5 * m3 * t2 +
                           YgaRs*XgAre * n1.c * w5 * m3 * t2 +
                           YgaRs*XgaRe * n1.c * w6 * m3 * t2 +
                           YgaRs*Xgare * n1.c * w6 * m3 * t2 +
                           YgaRs*YGARe * n1.c * w5 * m2 * t3 +
                           YgaRs*YGAre * n1.c * w5 * m2 * t3 +
                           YgaRs*YGaRe * n1.c * w6 * m2 * t3 +
                           YgaRs*YGare * n1.c * w6 * m2 * t3 +
                           YgaRs*YgARe * n1.c * w5 * m3 * t3 +
                           YgaRs*YgAre * n1.c * w5 * m3 * t3 +
                           YgaRs*YgaRe * n1.c * w6 * m3 * t3 +
                           YgaRs*Ygare * n1.c * w6 * m3 * t3 +
                           # row 16
                           Ygars*XGARe * n1.c * w5 * m2 * t2 +
                           Ygars*XGAre * n1.c * w5 * m2 * t2 +
                           Ygars*XGaRe * n1.c * w6 * m2 * t2 +
                           Ygars*XGare * n1.c * w6 * m2 * t2 +
                           Ygars*XgARe * n1.c * w5 * m3 * t2 +
                           Ygars*XgAre * n1.c * w5 * m3 * t2 +
                           Ygars*XgaRe * n1.c * w6 * m3 * t2 +
                           Ygars*Xgare * n1.c * w6 * m3 * t2 +
                           Ygars*YGARe * n1.c * w5 * m2 * t3 +
                           Ygars*YGAre * n1.c * w5 * m2 * t3 +
                           Ygars*YGaRe * n1.c * w6 * m2 * t3 +
                           Ygars*YGare * n1.c * w6 * m2 * t3 +
                           Ygars*YgARe * n1.c * w5 * m3 * t3 +
                           Ygars*YgAre * n1.c * w5 * m3 * t3 +
                           Ygars*YgaRe * n1.c * w6 * m3 * t3 +
                           Ygars*Ygare * n1.c * w6 * m3 * t3)
  
  if(return.fit == T){
    fits <- c(male.fit.bar, female.fit.bar)
    names(fits) <- c('male.fit','female.fit')
    stop(return(fits))
  }
  
  ##### XGARs.prime #####
  XGARs.prime <- (XGARs*XGARe * n2.c *w1*m1*t2 + 
                    XGARs*XGAre*.5*n2.c*w1*m1*t2 + 
                    XGARs*XGaRe *.5*n2.c*w2*m1*t2 +
                    XGARs*XGare * .5 * (r4.c*r5.c*r6.c +r4.c*r5*r6.c + 
                                          r4*r5.c*r6.c + r4*r5*r6.c)*n2.c*w2*m1*t2 + 
                    XGARs*XgARe * .5*n2.c*w1*m2*t2+
                    XGARs*XgAre*.5*(r4.c*r5.c*r6.c+r4.c*r5*r6+
                                      r4*r5.c*r6.c+r4*r5*r6)*n2.c*w1*m2*t2 +
                    XGARs*XgaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3+
                                      r1*r2.c*r3.c + r1*r2.c*r3) *
                    n2.c*w2*m2*t2 +
                    XGARs*Xgare *.5*(r4.c*r5.c*r6.c+r4*r5.c*r6.c)*
                    n2.c*w2*m2*t2 +
                    XGARs*YGARe * .5 *n1*w1*m1*t1 + 
                    XGARs*YGAre*.5*(r4.c*r5.c*r6.c+r4.c*r5*r6+r4*r5*r6.c + 
                                      r4*r5.c*r6)*n1 *w1*m1*t1 +
                    XGARs*YGaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3+
                                      r1*r2*r3.c+r1*r2*r3)*n1*w2*m1*t1+
                    XGARs*YGare*.5*(r4.c*r5.c*r6.c+r4*r5*r6.c)*n1*w2*m1*t1 +
                    XGARs*YgARe*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3.c + 
                                      r1.c*r2*r3+r1.c*r2.c*r3)*n1*w1*m2*t1 +
                    XGARs*YgAre*.5*(r4.c*r5.c*r6.c+r4.c*r5*r6)*n1*w1*m2*t1 +
                    XGARs*YgaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3)*
                    n1*w2*m2*t1 +
                    XGARs*Ygare*.5*r4.c*r5.c*r6.c*n1*w2*m2*t1 +
                    # row 2
                    XGArs*XGARe*.5*n2.c*w1*m1*t2 +
                    XGArs*XGaRe*.5*(r4.c*r5.c*r6+r4*r5*r6+ r4*r5.c*r6+
                                      r4.c*r5*r6)*n2.c*w2*m1*t2+
                    XGArs*XgARe*.5*(r4.c*r5.c*r6+r4.c*r5*r6.c + 
                                      r4*r5.c*r6+r4*r5*r6.c)*n2.c*w1*m2*t2 +
                    XGArs*XgaRe*.5*(r4.c*r5.c*r6+r4*r5.c*r6)*n2.c*w2*m2*t2 +
                    XGArs*YGARe*.5*(r4.c*r5.c*r6+r4*r5.c*r6.c+
                                      r4.c*r5*r6.c+r4*r5*r6)*n1*w1*m1*t1 +
                    XGArs*YGaRe*.5*(r4.c*r5.c*r6+r4*r5*r6)*n1*w2*m1*t1 +
                    XGArs*YgARe*.5*(r4.c*r5.c*r6+r4.c*r5*r6.c)*n1*w1*m2*t1 +
                    XGArs*YgaRe*.5*(r4.c*r5.c*r6)*n1*w2*m2*t1+
                    # row 3
                    XGaRs*XGARe*.5*n2.c*w2*m1*t2 +
                    XGaRs*XGAre*.5*(r4.c*r5*r6+r4*r5.c*r6+r4.c*r5.c*r6+
                                      r4*r5*r6)*n2.c*w2*m1*t2 +
                    XGaRs*XgARe*.5*(r1.c*r2*r3.c+r1.c*r2*r3+
                                      r1*r2*r3.c+r1*r2*r3)*n2.c*w2*m2*t2 +
                    XGaRs*XgAre*.5*(r4.c*r5*r6+r4*r5*r6)*n2.c*w2*m2*t2+
                    XGaRs*YGARe*.5*(r1*r2.c*r3+r1*r2.c*r3.c+
                                      r1.c*r2*r3.c + r1.c*r2*r3)*n1*w2*m1*t1 +
                    XGaRs*YGAre*.5*(r4.c*r5*r6+r4*r5.c*r6)*n1*w2*m1*t1 +
                    XGaRs*YgARe*.5*(r1.c*r2*r3+r1.c*r2*r3.c)*n1*w2*m2*t1 +
                    XGaRs*YgAre*.5*r4.c*r5*r6*n1*w2*m2*t1 +
                    # row 4
                    XGars*XGARe*.5*(r4.c*r5.c*r6.c+r4.c*r5*r6.c+ 
                                      r4*r5.c*r6.c+r4*r5*r6.c)*
                    n2.c*w2*m1*t2 +
                    XGars*XgARe*.5*(r4.c*r5*r6.c+r4*r5*r6.c)*n2.c*w2*m2*t2 +
                    XGars*YGARe*.5*(r4*r5.c*r6.c+r4.c*r5*r6.c)*n1*w2*m1*t1 +
                    XGars*YgARe*.5*r4.c*r5*r6.c*n1*w2*m2*t1 +
                    # row 5
                    XgARs*XGARe*.5*n2.c*w1*m2*t2 +
                    XgARs*XGAre*.5*(r4.c*r5*r6.c+r4.c*r5.c*r6+r4*r5.c*r6+
                                      r4*r5*r6.c)*n2.c*w1*m2*t2 +
                    XgARs*XGaRe*.5*(r1.c*r2*r3+r1.c*r2*r3.c+r1*r2*r3+
                                      r1*r2*r3.c)*n2.c*w2*m2*t2+
                    XgARs*XGare*.5*(r4.c*r5*r6.c+r4*r5*r6.c)*n2.c*w2*m2*t2 +
                    XgARs*YGARe*.5*(r1*r2*r3.c+r1*r2.c*r3.c+r1*r2*r3+
                                      r1*r2.c*r3)*n1*w1*m2*t1 +
                    XgARs*YGAre*.5*(r4*r5*r6.c+r4*r5.c*r6)*n1*w1*m2*t1 +
                    XgARs*YGaRe*.5*(r1*r2*r3+r1*r2*r3.c)*n1*w2*m2*t1 +
                    XgARs*YGare*.5*r4*r5*r6.c*n1*w2*m2*t1 +
                    # row 6
                    XgArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + 
                                          r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XgArs*XGaRe * .5 * (r4*r5*r6 +r4.c*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgArs*YGARe * .5 * (r4*r5*r6 + r4*r5.c*r6.c) *n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + 
                                          r1.c*r2.c*r3 + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGAre * .5 * (r4*r5.c*r6 + r4.c*r5.c*r6) * n2.c * w2 *m2 *t2 +
                    XgaRs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 *w2 *m2 *t1 + 
                    XgaRs*YGAre * .5 * r4*r5.c*r6 * n1 *w2 * m2 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2+
                    Xgars*YGARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 *t1 + 
                    # row 9
                    YGARs*XGARe * .5 * n1 * w1 * m1 * t1 +
                    YGARs*XGAre * .5 * (r4*r5.c*r6.c + r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 + 
                    YGARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m1 * t1 + 
                    YGARs*XGare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 + 
                    YGARs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3 +r1*r2*r3.c) * n1 * w1 * m2 * t1 +
                    YGARs*XgAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YGARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6 + r4*r5*r6.c + r4.c*r5*r6) * n1 * w1 * m1 * t1 +
                    YGArs*XGaRe * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m1 * t1 +
                    YGArs*XgARe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * (r4*r5.c*r6) * n1 * w2 * m2 * t1 + 
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    YGaRs*XGAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGaRs*XgARe * .5 * (r1*r2*r3 + r1*r2*r3.c) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGars*XgARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 13
                    YgARs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3 + r1.c*r2*r3.c) * n1 *w1 * m2 * t1 +
                    YgARs*XGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 + 
                    YgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 14
                    YgArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 + 
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1) / male.fit.bar
  
  ##### XGArs.prime #####
  XGArs.prime <- (XGARs*XGAre * .5 * n2.c * w1 * m1 * t2 +
                    XGARs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5*r6 + r4*r5.c*r6 + r4.c*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGARs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n2.c * w1 * m2 * t2 +
                    XGARs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    XGARs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGARs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 + 
                    XGARs*Ygare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XGARe * .5 * n2.c * w1 * m1 * t2 +
                    XGArs*XGAre * n2.c * w1 * m1 * t2 +
                    XGArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6.c) * n2.c * w2 * m1 * t2 +
                    XGArs*XGare * .5 * n2.c * w2 * m1 * t2 +
                    XGArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XGArs*XgAre * .5 * n2.c * w1 * m2 * t2 +
                    XGArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XGArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    XGArs*YGAre * .5 * n1 * w1 * m1 * t1 +
                    XGArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    XGArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XGArs*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8*r9 + r7.c*r8*r9.c + r7.c*r8.c*r9) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6.c) * n2.c * w2 * m1 * t2 +
                    XGaRs*XgAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGaRs*YGAre * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    XGaRs*YgAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*XGARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6 + r4.c*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGars*XGAre * .5 * n2.c * w2 * m1 * t2 +
                    XGars*XgARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2+
                    XGars*XgAre * .5 * (r7.c*r8*r9 + r7.c*r8*r9.c + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XGars*YGARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 + 
                    XGars*YGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    XGars*YgARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XgARs*XGare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgARs*YGAre * .5 * (r4*r5.c*r6.c +r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XgARs*YGare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n2.c * w1 * m2 * t2 +
                    XgArs*XGAre * .5 * n2.c * w1 * m2 * t2 +
                    XgArs*XGaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XgArs*XGare * .5 * (r7.c*r8*r9 + r7.c*r8*r9.c + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XgArs*YGARe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 + 
                    XgArs*YGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgaRs*YGAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    Xgars*XGAre * .5 * (r7.c*r8.c*r9.c + r7*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9) * n2.c * w2 * m2 * t2 + 
                    Xgars*YGARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 9
                    YGARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    YGARs*XGare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 + 
                    YGARs*XgAre * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 + 
                    YGARs*Xgare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    YGArs*XGAre * .5 * n1 * w1 * m1 * t1 +
                    YGArs*XGaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    YGArs*XGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    YGArs*XgARe * .5 * (r4*r5.c*r6.c +r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YGArs*XgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGaRs*XgAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGars*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    YGars*XgARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 13
                    YgARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 14
                    YgArs*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 + 
                    YgArs*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8*r9 + r7.c*r8*r9.c + r7.c*r8.c*r9) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 15
                    YgaRs*XGAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1) / male.fit.bar
  
  ##### XGaRs.prime #####
  XGaRs.prime <- (XGARs*XGaRe * .5 * n2.c * w2 * m1 * t2 +
                    XGARs*XGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XGARs*Xgare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XGARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m1 * t1 +
                    XGARs*YGare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 + 
                    XGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 + 
                    XGARs*Ygare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 + 
                    # row 2
                    XGArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m1 * t2 +
                    XGArs*XgaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGArs*YGaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    XGArs*YgaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XGARe * .5 * n2.c * w2 * m1 * t2 +
                    XGaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m1 * t2 +
                    XGaRs*XGaRe * n2.c * w3 * m1 * t2 +
                    XGaRs*XGare * .5 * n2.c * w3 * m1 * t2 +
                    XGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XGaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGaRs*XgaRe * .5 * n2.c * w3 * m2 * t2 +
                    XGaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XGaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    XGaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGaRs*YGaRe * .5 * n1 * w3 * m1 * t1 +
                    XGaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m1 * t1 +
                    XGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XGaRs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGars*XGaRe * .5 * n2.c * w3 * m1 * t2 +
                    XGars*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGars*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n2.c * w3 * m2 * t2 +
                    XGars*YGARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGars*YGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    XGars*YgARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 5
                    XgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XgARs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XgARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgArs*YGaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGaRe * .5 * n2.c * w3 * m2 * t2 +
                    XgaRs*XGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n2.c * w3 * m2 * t2   +
                    XgaRs*YGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    Xgars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    Xgars*YGARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 9
                    YGARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    YGARs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGARs*XgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGArs*XgaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m1 * t1 +
                    YGaRs*XGAre * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    YGaRs*XGaRe * .5 * n1 * w3 * m1 * t1 +
                    YGaRs*XGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    YGaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 + 
                    YGars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m1 * t1 +
                    YGars*XgARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 14
                    YgArs*XGaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 + 
                    YgaRs*XGAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 + 
                    Ygars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1) / male.fit.bar
  
  ##### XGars.prime #####
  XGars.prime <- (XGARs*XGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + 
                                        r4*r5*r6.c) * n2.c * w2 * m1 * t2 + 
                    XGARs*Xgare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGARs*YGare * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    XGARs*Ygare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGArs*XGare * .5 * n2.c * w2 * m1 * t2 +
                    XGArs*XgaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XGArs*Xgare * .5 * (r7.c*r8*r9.c + r7*r8*r9.c + r7.c*r8*r9 + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XGArs*YGaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    XGArs*YGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m1 * t1 + 
                    XGArs*YgaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m1 * t2 +
                    XGaRs*XGare * .5 * n2.c * w3 * m1 * t2 +
                    XGaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGaRs*Xgare * .5 * (r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6.c + r4*r5.c*r6) * n2.c * w3 * m2 * t2 +
                    XGaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGaRs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    XGaRs*YgAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*XGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6.c) * n2.c * w2 * m1 * t2 + 
                    XGars*XGAre * .5 * n2.c * w2 * m1 * t2 +
                    XGars*XGaRe * .5 * n2.c * w3 * m1 * t2 +
                    XGars*XGare * n2.c * w3 * m1 * t2 +
                    XGars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XGars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XGars*Xgare * .5 * n2.c * w3 * m2 * t2 +
                    XGars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 + 
                    XGars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 + 
                    XGars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m1 * t1 +
                    XGars*YGare * .5 * n1 * w3 * m1 * t1 +
                    XGars*YgARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 + 
                    XGars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    XGars*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 5 
                    XgARs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgARs*YGare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 + 
                    # row 6
                    XgArs*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XgArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XgArs*YGaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 7 
                    XgaRs*XGAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XgaRs*YGAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    Xgars*XGAre * .5 * (r7.c*r8*r9.c + r7*r8*r9.c + r7.c*r8*r9 + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    Xgars*XGaRe * .5 * (r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6.c + r4*r5.c*r6) * n2.c * w3 * m2 * t2 +
                    Xgars*XGare * .5 * n2.c * w3 * m2 * t2 +
                    Xgars*YGARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    Xgars*YGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 9
                    YGARs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 + 
                    YGARs*Xgare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 10 
                    YGArs*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 + 
                    YGArs*XgaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XGAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    YGaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m1 * t1 +
                    YGaRs*XgAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    YGars*XGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m1 * t1 + 
                    YGars*XGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    YGars*XGare * .5 * n1 * w3 * m1 * t1 +
                    YGars*XgARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 + 
                    YGars*XgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    YGars*Xgare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 13
                    YgARs*XGare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 + 
                    # row 14
                    YgArs*XGaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 15
                    YgaRs*XGAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    Ygars*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1) / 
    male.fit.bar
  
  ##### XgARs.prime #####
  XgARs.prime <- (XGARs*XgARe * .5 * n2.c * w1 * m2 * t2 + 
                    XGARs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c + r4*r5.c*r6) * n2.c * w1 * m2 * t2 +
                    XGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XGARs*Xgare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGARs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    XGARs*YgAre * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XGARs*YgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 + 
                    XGARs*Ygare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XGArs*XgaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XGArs*YgARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XGaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGars*YgARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*XGARe * .5 * n2.c * w1 * m2 * t2 + 
                    XgARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XgARs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgARs*XgARe * n2.c * w1 * m3 * t2 +
                    XgARs*XgAre * .5 * n2.c * w1 * m3 * t2 +
                    XgARs*XgaRe * .5 * n2.c * w2 * m3 * t2 +
                    XgARs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgARs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 +  r1.c*r2*r3) * n1 * w1 * m2 * t1 +
                    XgARs*YGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgARs*YgARe * .5 * n1 * w1 * m3 * t1 +
                    XgARs*YgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m3 * t1 +
                    XgARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*XGARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w1 * m2 * t2 +
                    XgArs*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6)* n2.c * w2 * m2 * t2 +
                    XgArs*XgARe * .5 * n2.c * w1 * m3 * t2 +
                    XgArs*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgArs*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgArs*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgaRs*XgARe * .5 * n2.c * w2 * m3 * t2 +
                    XgaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XgaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m3 * t1 +
                    XgaRs*YgAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    Xgars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    Xgars*YGARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 + 
                    Xgars*YgARe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    # row 9 
                    YGARs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 +  r1.c*r2*r3) * n1 * w1 * m2 * t1 +
                    YGARs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 + 
                    # row 10
                    YGArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 12
                    YGars*XgARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 13
                    YgARs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    YgARs*XGare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YgARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgARs*XgARe * .5 * n1 * w1 * m3 * t1 +
                    YgARs*XgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    YgARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    # row 14
                    YgArs*XGARe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m3 * t1 +
                    YgArs*XgaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 + 
                    YgaRs*XGAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1) / male.fit.bar
  
  ##### XgArs.prime #####
  XgArs.prime <- (XGARs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XGARs*Xgare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 *  m2 * t2 +
                    XGARs*YgAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XGARs*Ygare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w1 * m2 * t2 +
                    XGArs*XgAre * .5 * n2.c * w1 * m2 * t2 +
                    XGArs*XgaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGArs*Xgare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XGArs*YgARe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XGArs*YgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGaRs*YgAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XGars*YgARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*XGAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w1 * m2 * t2 +
                    XgARs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XgARs*XgAre * .5 * n2.c * w1 * m3 * t2 +
                    XgARs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgARs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w1 * m2 * t2 +
                    XgArs*XGAre * .5 * n2.c * w1 * m2 * t2 +
                    XgArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XgArs*XgARe * .5 * n2.c * w1 * m3 * t2 +
                    XgArs*XgAre * n2.c * w1 * m3 * t2 +
                    XgArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgArs*Xgare * .5 * n2.c * w2 * m3 * t2 +
                    XgArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XgArs*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    XgArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    XgArs*YgAre * .5 * n1 * w1 * m3 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    XgArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*XGAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XgaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgaRs*YGAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YgAre * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 *  m2 * t2 +
                    Xgars*XGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    Xgars*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    Xgars*XgAre * .5 * n2.c * w2 * m3 * t2 +
                    Xgars*YGARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m3 * t1 +
                    Xgars*YgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m3 * t1 +
                    # row 9
                    YGARs*XgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGArs*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XgAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 12
                    YGars*XgARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 13
                    YgARs*XGAre * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgARs*XgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 14
                    YgArs*XGARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YgArs*XGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YgArs*XgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    YgArs*XgAre * .5 * n1 * w1 * m3 * t1 +
                    YgArs*XgaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    YgArs*Xgare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m3 * t1 +
                    # row 15
                    YgaRs*XGAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    Ygars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1) / male.fit.bar
  
  ##### XgaRs.prime #####
  XgaRs.prime <- (XGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + 
                                        r1*r2.c*r3.c + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XGARs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGARs*Ygare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGArs*YgaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XGaRs*XgAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGaRs*XgaRe * .5 * n2.c * w3 * m2 * t2 +
                    XGaRs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w3 * m2 * t2 +
                    XGaRs*YgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XGaRs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*XgARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XGars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XGars*YgARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 5
                    XgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2.c * w2 * m2 * t2 +
                    XgARs*XGare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgARs*XgaRe * .5 * n2.c * w2 * m3 * t2 +
                    XgARs*Xgare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XgARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*XGaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XgArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgArs*YGaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YgaRe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGaRe * .5 * n2.c * w3 * m2 * t2 +
                    XgaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XgaRs*XgARe * .5 * n2.c * w2 * m3 * t2 +
                    XgaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgaRs*XgaRe * n2.c * w3 * m3 * t2 + 
                    XgaRs*Xgare * .5 * n2.c * w3 * m3 * t2 + 
                    XgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    XgaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    XgaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    XgaRs*YgaRe * .5 * n1 * w3 * m3 * t1 +
                    XgaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m3 * t1 +
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    Xgars*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w3 * m2 * t2 +
                    Xgars*XgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    Xgars*XgaRe * .5 * n2.c * w3 * m3 * t2 + 
                    Xgars*YGARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    Xgars*YgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m3 * t1 +
                    # row 9
                    YGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XgaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 12
                    YGars*XgARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 14
                    YgArs*XGaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    YgaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m3 * t1 +
                    YgaRs*XgAre * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgaRs*XgaRe * .5 * n1 * w3 * m3 * t1 +
                    YgaRs*Xgare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m3 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    Ygars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m3 * t1) / male.fit.bar
  
  ##### Xgars.prime #####
  Xgars.prime <- (XGARs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    XGARs*Ygare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*XgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XGArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    XGArs*YgaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*XgAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XGaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    XGaRs*YgAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*XgARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XGars*XgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XGars*XgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w3 * m2 * t2 +
                    XGars*Xgare * .5 * n2.c * w3 * m2 * t2 +
                    XGars*YgARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 + 
                    XGars*YgAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 + 
                    XGars*YgaRe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 + 
                    XGars*Ygare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 + 
                    # row 5
                    XgARs*XGare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m2 * t2 +
                    XgARs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    XgARs*YGare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 + 
                    XgARs*Ygare * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 + 
                    # row 6
                    XgArs*XGaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n2.c * w2 * m2 * t2 +
                    XgArs*XGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n2.c * w2 * m2 * t2 +
                    XgArs*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgArs*Xgare * .5 * n2.c * w2 * m3 * t2 +
                    XgArs*YGaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 + 
                    XgArs*YGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 + 
                    XgArs*YgaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 + 
                    XgArs*Ygare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m3 * t1 + 
                    # row 7
                    XgaRs*XGAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n2.c * w2 * m2 * t2 +
                    XgaRs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w3 * m2 * t2 +
                    XgaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n2.c * w2 * m3 * t2 +
                    XgaRs*Xgare * .5 * n2.c * w3 * m3 * t2 +
                    XgaRs*YGAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 + 
                    XgaRs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 + 
                    XgaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 + 
                    XgaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m3 * t1 + 
                    # row 8
                    Xgars*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n2.c * w2 * m2 * t2 +
                    Xgars*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n2.c * w2 * m2 * t2 +
                    Xgars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n2.c * w3 * m2 * t2 +
                    Xgars*XGare * .5 * n2.c * w3 * m2 * t2 +
                    Xgars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n2.c * w2 * m3 * t2 +
                    Xgars*XgAre * .5 * n2.c * w2 * m3 * t2 +
                    Xgars*XgaRe * .5 * n2.c * w3 * m3 * t2 +
                    Xgars*Xgare * n2.c * w3 * m3 * t2 +
                    Xgars*YGARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    Xgars*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    Xgars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    Xgars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m3 * t1 +
                    Xgars*Ygare * .5 * n1 * w3 * m3 * t1 +
                    # row 9
                    YGARs*Xgare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 10
                    YGArs*XgaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 + 
                    YGArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 11
                    YGaRs*XgAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 + 
                    YGaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 12
                    YGars*XgARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 + 
                    YGars*XgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 + 
                    YGars*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 + 
                    YGars*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 13
                    YgARs*XGare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 + 
                    YgARs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 14
                    YgArs*XGaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 + 
                    YgArs*XgaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 + 
                    YgArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    # row 15
                    YgaRs*XGAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 + 
                    YgaRs*XgAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 + 
                    YgaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w3 * m3 * t1 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    Ygars*XGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 + 
                    Ygars*XgARe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 + 
                    Ygars*XgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m3 * t1 + 
                    Ygars*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m3 * t1 + 
                    Ygars*Xgare * .5 * n1 * w3 * m3 * t1) / male.fit.bar
  
  ##### YGARs.prime #####
  YGARs.prime <- (XGARs*YGARe * .5 * n1 * w1 * m1 * t1 +
                    XGARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    XGARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m1 * t1 +
                    XGARs*YGare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGARs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    XGARs*YgAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XGARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGARs*Ygare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    XGArs*YGaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    XGArs*YgARe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    XGaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGaRs*YgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGars*YgARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 5 
                    XgARs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w1 * m2 * t1 +
                    XgARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 8 
                    Xgars*YGARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 9
                    YGARs*XGARe * .5 * n1 * w1 * m1 * t1 +
                    YGARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    YGARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    YGARs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGARs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w1 * m2 * t1 +
                    YGARs*XgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGARs*YGARe * n1 * w1 * m1 * t2 +
                    YGARs*YGAre * .5 * n1 * w1 * m1 * t2 +
                    YGARs*YGaRe * .5 * n1 * w2 * m1 * t2 +
                    YGARs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGARs*YgARe * .5 * n1 * w1 * m2 * t2 +
                    YGARs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    YGArs*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGArs*XgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGArs*YGARe * .5 * n1 * w1 * m1 * t2 +
                    YGArs*YGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGArs*YgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m1 * t1 +
                    YGaRs*XGAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    YGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YGaRs*YGARe * .5 * n1 * w2 * m1 * t2 +
                    YGaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGars*XgARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGars*YgARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 13
                    YgARs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    YgARs*XGAre * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 +
                    YgARs*XGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgARs*YGARe * .5 * n1 * w1 * m2 * t2 +
                    YgARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YgaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2) / male.fit.bar
  
  ##### YGArs.prime #####
  YGArs.prime <- (XGARs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    XGARs*YGare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    XGARs*YgAre * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 +
                    XGARs*Ygare * .5 * (r4*r5.c*r6) * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    XGArs*YGAre * .5 * n1 * w1 * m1 * t1 +
                    XGArs*YGaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    XGArs*YGare * .5 * (r7.c*r8*r9.c + r7*r8.c*r9.c + r7.c*r8*r9 + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    XGArs*YgARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XGArs*YgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGaRs*YgAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*YGARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    XGars*YgARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*YGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgArs*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*YGAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 9
                    YGARs*XGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m1 * t1 +
                    YGARs*XGare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGARs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGARs*YGAre * .5 * n1 * w1 * m1 * t2 +
                    YGARs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGARs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6 + r4*r5.c*r6) * n1 * w1 * m1 * t1 +
                    YGArs*XGAre * .5 * n1 * w1 * m1 * t1 +
                    YGArs*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    YGArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YGArs*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YGArs*YGARe * .5 * n1 * w1 * m1 * t2 +
                    YGArs*YGAre * n1 * w1 * m1 * t2 +
                    YGArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGArs*YGare * .5 * n1 * w2 * m1 * t2 +
                    YGArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YGArs*YgAre * .5 * n1 * w1 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YGArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XGAre * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    YGars*XGAre * .5 * (r7.c*r8*r9.c + r7*r8.c*r9.c + r7.c*r8*r9 + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    YGars*XgARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    YGars*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGars*YGAre * .5 * n1 * w2 * m1 * t2 +
                    YGars*YgARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YGars*YgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    # row 13
                    YgARs*XGAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgARs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w1 * m2 * t1 +
                    YgArs*XGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    YgArs*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YgArs*YGAre * .5 * n1 * w1 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YgArs*YGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * (r4*r5.c*r6) * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    Ygars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2) / 
    male.fit.bar
  
  ##### YGaRs.prime #####
  YGaRs.prime <- (XGARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    XGARs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    XGARs*YgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XGARs*Ygare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGArs*YgaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m1 * t1 +
                    XGaRs*YGAre * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGaRs*YGaRe * .5 * n1 * w3 * m1 * t1 +
                    XGaRs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    XGaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XGaRs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*YGARe * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m1 * t1 +
                    XGars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6 + r4.c*r5*r6) * n1 * w3 * m1 * t1 +
                    XGars*YgARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 5
                    XgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*YGaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 8 
                    Xgars*YGARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 9
                    YGARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m1 * t1 +
                    YGARs*XGare * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m1 * t1 +
                    YGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YGARs*YGaRe * .5 * n1 * w2 * m1 * t2 +
                    YGARs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XGaRe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t2 +
                    YGArs*YgaRe * .5 * (r4*r5*r6.c + r4.c*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m1 * t1 +
                    YGaRs*XGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGaRs*XGaRe * .5 * n1 * w3 * m1 * t1 +
                    YGaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6 + r4.c*r5*r6) * n1 * w3 * m1 * t1 +
                    YGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    YGaRs*YGARe * .5 * n1 * w2 * m1 * t2 +
                    YGaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m1 * t2 +
                    YGaRs*YGaRe * n1 * w3 * m1 * t2 +
                    YGaRs*YGare * .5 * n1 * w3 * m1 * t2 +
                    YGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YGaRs*YgaRe * .5 * n1 * w3 * m2 * t2 +
                    YGaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6)* n1 * w3 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGars*XGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    YGars*XgARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    YGars*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGars*YGaRe * .5 * n1 * w3 * m1 * t2 +
                    YGars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YGars*YgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    YgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YgaRs*YGAre * .5 * (r4*r5*r6.c + r4.c*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YgaRs*YGaRe * .5 * n1 * w3 * m2 * t2 +
                    YgaRs*YGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    Ygars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w3 * m2 * t2) / male.fit.bar
  
  ##### YGars.prime #####
  YGars.prime <- (XGARs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGARs*Ygare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YGaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 + # 2nd r5 was r5.c
                    XGArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    XGArs*YgaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YGAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    XGaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6 + r4.c*r5*r6) * n1 * w3 * m1 * t1 +
                    XGaRs*YgAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*YGARe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    XGars*YGAre * .5 * (r7*r8.c*r9.c + r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    XGars*YGaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    XGars*YGare * .5 * n1 * w3 * m1 * t1 +
                    XGars*YgARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 + 
                    XGars*Ygare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 + 
                    # row 5
                    XgARs*YGare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 6
                    XgArs*YGaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 7
                    XgaRs*YGAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    Xgars*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 9
                    YGARs*XGare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGARs*Xgare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGARs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGARs*Ygare * .5 * (r4*r5*r6.c + r4.c*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XGaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m1 * t1 +
                    YGArs*XGare * .5 * (r7*r8.c*r9.c + r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8.c*r9) * n1 * w2 * m1 * t1 +
                    YGArs*XgaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    YGArs*YGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGArs*YGare * .5 * n1 * w2 * m1 * t2 +
                    YGArs*YgaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YGArs*Ygare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XGAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t1 +
                    YGaRs*XGare * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w3 * m1 * t1 +
                    YGaRs*XgAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    YGaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m1 * t2 +
                    YGaRs*YGare * .5 * n1 * w3 * m1 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YGaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m1 * t1 +
                    YGars*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m1 * t1 +
                    YGars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6 + r4.c*r5*r6) * n1 * w3 * m1 * t1 +
                    YGars*XGare * .5 * n1 * w3 * m1 * t1 +
                    YGars*XgARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    YGars*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    YGars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m1 * t2 +
                    YGars*YGAre * .5 * n1 * w2 * m1 * t2 +
                    YGars*YGaRe * .5 * n1 * w3 * m1 * t2 +
                    YGars*YGare * n1 * w3 * m1 * t2 +
                    YGars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YGars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    YGars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    YGars*Ygare * .5 * n1 * w3 * m2 * t2 +
                    # row 13
                    YgARs*XGare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgARs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YgArs*YGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YgArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 + 
                    YgaRs*YGAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YgaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    Ygars*XGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 + 
                    Ygars*YGARe * .5 * (r4*r5*r6.c + r4.c*r5*r6.c) * n1 * w2 * m2 * t2 +
                    Ygars*YGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    Ygars*YGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    Ygars*YGare * .5 * n1 * w3 * m2 * t2) / male.fit.bar
  
  ##### YgARs.prime #####
  YgARs.prime <- (XGARs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * 
                    n1 * w1 * m2 * t1 +
                    XGARs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XGARs*Ygare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*YgARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    XgARs*YGAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XgARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgARs*YgARe * .5 * n1 * w1 * m3 * t1 +
                    XgARs*YgAre * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    XgARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*YGARe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XgaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    XgaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 9
                    YGARs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w1 * m2 * t1 +
                    YGARs*XgAre * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGARs*XgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGARs*YgARe * .5 * n1 * w1 * m2 * t2 +
                    YGARs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XgARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGArs*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 13
                    YgARs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w1 * m2 * t1 +
                    YgARs*XGAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgARs*XgARe * .5 * n1 * w1 * m3 * t1 +
                    YgARs*XgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    YgARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgARs*YGARe * .5 * n1 * w1 * m2 * t2 +
                    YgARs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YgARs*YgARe * n1 * w1 * m3 * t2 +
                    YgARs*YgAre * .5 * n1 * w1 * m3 * t2 +
                    YgARs*YgaRe * .5 * n1 * w2 * m3 * t2 +
                    YgARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XgARe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    YgArs*XgaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    YgArs*YGARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YgArs*YgARe * .5 * n1 * w1 * m3 * t2 +
                    YgArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m3 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    YgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t2 +
                    YgaRs*YGAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YgaRs*YgARe * .5 * n1 * w2 * m3 * t2 +
                    YgaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    Ygars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t2) / 
    male.fit.bar
  
  
  ##### YgArs.prime #####
  YgArs.prime <- (XGARs*YgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    XGARs*Ygare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XGArs*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8*r9 + r7.c*r8.c*r9 + r7.c*r8*r9.c) * n1 * w1 * m2 * t1 +
                    XGArs*YgaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YgAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 4
                    XGars*YgARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 5
                    XgARs*YGAre * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    XgARs*YGare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgARs*YgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*YGARe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    XgArs*YGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    XgArs*YGaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    XgArs*YgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    XgArs*YgAre * .5 * n1 * w1 * m3 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    XgArs*Ygare * .5 * (r7.c*r8*r9.c + r7*r8.c*r9.c + r7*r8.c*r9 + r7.c*r8*r9) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*YGAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    Xgars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    # row 9
                    YGARs*XgAre * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t1 +
                    YGARs*Xgare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGARs*YgAre * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XgARe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YGArs*XgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w1 * m2 * t1 +
                    YGArs*XgaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    YGArs*YgARe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YGArs*YgAre * .5 * n1 * w1 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YGArs*Ygare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XgAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YGars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YGars*YgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    # row 13
                    YgARs*XGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w1 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgARs*XgAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    YgARs*YGAre * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YgARs*YgAre * .5 * n1 * w1 * m3 * t2 +
                    YgARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w1 * m2 * t1 +
                    YgArs*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8*r9 + r7.c*r8.c*r9 + r7.c*r8*r9.c) * n1 * w1 * m2 * t1 +
                    YgArs*XGaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YgArs*XgARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w1 * m3 * t1 +
                    YgArs*XgAre * .5 * n1 * w1 * m3 * t1 +
                    YgArs*XgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgArs*Xgare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    YgArs*YGARe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n1 * w1 * m2 * t2 +
                    YgArs*YGAre * .5 * n1 * w1 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YgArs*YGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    YgArs*YgARe * .5 * n1 * w1 * m3 * t2 +
                    YgArs*YgAre * n1 * w1 * m3 * t2 +
                    YgArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c +  r4*r5*r6.c) * n1 * w2 * m3 * t2 +
                    YgArs*Ygare * .5 * n1 * w2 * m3 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5*r6.c + r4*r5.c*r6.c) * n1 * w2 * m3 * t1 +
                    YgaRs*YGAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YgaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c +  r4*r5*r6.c) * n1 * w2 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    Ygars*XgAre * .5 * (r7.c*r8*r9.c + r7*r8.c*r9.c + r7*r8.c*r9 + r7.c*r8*r9) * n1 * w2 * m3 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    Ygars*YGAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    Ygars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    Ygars*YgAre * .5 * n1 * w2 * m3 * t2) / male.fit.bar
  
  ##### YgaRs.prime #####
  YgaRs.prime <- (XGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XGARs*Ygare * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YgaRe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    XGaRs*YgAre * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XGaRs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*YgARe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 5
                    XgARs*YGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    XgARs*YGare * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XgARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*YGaRe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    XgaRs*YGAre * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    XgaRs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    XgaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m3 * t1 +
                    XgaRs*YgAre * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    XgaRs*YgaRe * .5 * n1 * w3 * m3 * t1 +
                    XgaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w3 * m3 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    Xgars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m3 * t1 +
                    # row 9
                    YGARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YGARs*Xgare * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YGARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XgaRe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m2 * t1 +
                    YGaRs*XgAre * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGaRs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w3 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m2 * t1 +
                    YGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3)  * n1 * w2 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YGaRs*YgaRe * .5 * n1 * w3 * m2 * t2 +
                    YGaRs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    YGars*YgARe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YGars*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m2 * t1 +
                    YgARs*XGare * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YgARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w2 * m3 * t1 +
                    YgARs*Xgare * .5 * (r4.c*r5*r6 + r4*r5.c*r6) * n1 * w2 * m3 * t1 +
                    YgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3)  * n1 * w2 * m2 * t2 +
                    YgARs*YGare * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YgARs*YgaRe * .5 * n1 * w2 * m3 * t2 +
                    YgARs*Ygare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgArs*XgaRe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgArs*YGaRe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YgArs*YgaRe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1 * w2 * m2 * t1 +
                    YgaRs*XGAre * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YgaRs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1 * w3 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    YgaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1 * w2 * m3 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgaRs*XgaRe * .5 * n1 * w3 * m3 * t1 +
                    YgaRs*Xgare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6.c + r4*r5.c*r6) * n1 * w3 * m3 * t1 +
                    YgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1 * w2 * m2 * t2 +
                    YgaRs*YGAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    YgaRs*YGaRe * .5 * n1 * w3 * m2 * t2 +
                    YgaRs*YGare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6 + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    YgaRs*YgARe * .5 * n1 * w2 * m3 * t2 +
                    YgaRs*YgAre * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t2 +
                    YgaRs*YgaRe * n1 * w3 * m3 * t2 +
                    YgaRs*Ygare * .5 * n1 * w3 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    Ygars*XgaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6 + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w3 * m3 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    Ygars*YGaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    Ygars*YgARe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    Ygars*YgaRe * .5 * n1 * w3 * m3 * t2) / male.fit.bar
  
  ##### Ygars.prime #####
  Ygars.prime <- (XGARs*Ygare * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    # row 2
                    XGArs*YgaRe * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XGArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    # row 3
                    XGaRs*YgAre * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    XGaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    # row 4
                    XGars*YgARe * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    XGars*YgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    XGars*YgaRe * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    XGars*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    # row 5
                    XgARs*YGare * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    XgARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    # row 6
                    XgArs*YGaRe * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    XgArs*YGare * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    XgArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    XgArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    # row 7 
                    XgaRs*YGAre * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    XgaRs*YGare * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    XgaRs*YgAre * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m3 * t1 +
                    XgaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m3 * t1 +
                    # row 8
                    Xgars*YGARe * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    Xgars*YGAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    Xgars*YGaRe * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    Xgars*YGare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 +
                    Xgars*YgARe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    Xgars*YgAre * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m3 * t1 +
                    Xgars*YgaRe * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w3 * m3 * t1 +
                    Xgars*Ygare * .5 * n1 * w3 * m3 * t1 +
                    # row 9
                    YGARs*Xgare * .5 * r4*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    YGARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    # row 10
                    YGArs*XgaRe * .5 * r4*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YGArs*Xgare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t1 +
                    YGArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YGArs*Ygare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    # row 11
                    YGaRs*XgAre * .5 * r4*r5*r6 * n1 * w2 * m2 * t1 +
                    YGaRs*Xgare * .5 * (r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t1 +
                    YGaRs*YgAre * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YGaRs*Ygare * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r4*r5*r6.c * n1 * w2 * m2 * t1 +
                    YGars*XgAre * .5 * (r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t1 +
                    YGars*XgaRe * .5 * (r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m2 * t1 +
                    YGars*Xgare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w3 * m2 * t1 +
                    YGars*YgARe * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YGars*YgAre * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    YGars*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    YGars*Ygare * .5 * n1 * w3 * m2 * t2 +
                    # row 13
                    YgARs*XGare * .5 * r4.c*r5*r6.c * n1 * w2 * m2 * t1 +
                    YgARs*Xgare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t1 +
                    YgARs*YGare * .5 * (r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w2 * m2 * t2 +
                    YgARs*Ygare * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r4.c*r5*r6 * n1 * w2 * m2 * t1 +
                    YgArs*XGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m2 * t1 +
                    YgArs*XgaRe * .5 * (r4*r5.c*r6 + r4.c*r5*r6) * n1 * w2 * m3 * t1 +
                    YgArs*Xgare * .5 * (r7*r8.c*r9.c + r7*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w2 * m3 * t1 +
                    YgArs*YGaRe * .5 * (r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m2 * t2 +
                    YgArs*YGare * .5 * (r7.c*r8*r9.c + r7.c*r8*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m2 * t2 +
                    YgArs*YgaRe * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    YgArs*Ygare * .5 * n1 * w2 * m3 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r4.c*r5.c*r6 * n1 * w2 * m2 * t1 +
                    YgaRs*XGare * .5 * (r4.c*r5.c*r6 + r4.c*r5*r6.c) * n1 * w3 * m2 * t1 +
                    YgaRs*XgAre * .5 * (r4.c*r5.c*r6 + r4*r5*r6) * n1 * w2 * m3 * t1 +
                    YgaRs*Xgare * .5 * (r4*r5.c*r6.c + r4.c*r5*r6.c + r4.c*r5.c*r6 + r4*r5*r6) * n1 * w3 * m3 * t1 +
                    YgaRs*YGAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6) * n1 * w2 * m2 * t2 +
                    YgaRs*YGare * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6.c + r4*r5*r6.c) * n1 * w3 * m2 * t2 +
                    YgaRs*YgAre * .5 * (r4.c*r5.c*r6 + r4*r5.c*r6 + r4.c*r5*r6 + r4*r5*r6) * n1 * w2 * m3 * t2 +
                    YgaRs*Ygare * .5 * n1 * w3 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r4.c*r5.c*r6.c * n1 * w2 * m2 * t1 +
                    Ygars*XGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9) * n1 * w2 * m2 * t1 +
                    Ygars*XGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6) * n1 * w3 * m2 * t1 +
                    Ygars*XGare * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7.c*r8*r9.c + r7.c*r8*r9) * n1 * w3 * m2 * t1 +
                    Ygars*XgARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c) * n1 * w2 * m3 * t1 +
                    Ygars*XgAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8*r9.c + r7*r8*r9) * n1 * w2 * m3 * t1 +
                    Ygars*XgaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6 + r4*r5*r6.c) * n1 * w3 * m3 * t1 +
                    Ygars*Xgare * .5 * n1 * w3 * m3 * t1 +
                    Ygars*YGARe * .5 * (r4.c*r5.c*r6.c + r4*r5.c*r6.c) * n1 * w2 * m2 * t2 +
                    Ygars*YGAre * .5 * (r7.c*r8.c*r9.c + r7.c*r8.c*r9 + r7*r8.c*r9.c + r7*r8.c*r9) * n1 * w2 * m2 * t2 +
                    Ygars*YGaRe * .5 * (r4.c*r5.c*r6.c + r4.c*r5*r6 + r4*r5.c*r6.c + r4*r5*r6) * n1 * w3 * m2 * t2 +
                    Ygars*YGare * .5 * n1 * w3 * m2 * t2 +
                    Ygars*YgARe * .5 * (r4.c*r5.c*r6.c + r4*r5*r6.c + r4*r5.c*r6.c + r4.c*r5*r6.c) * n1 * w2 * m3 * t2 +
                    Ygars*YgAre * .5 * n1 * w2 * m3 * t2 +
                    Ygars*YgaRe * .5 * n1 * w3 * m3 * t2 +
                    Ygars*Ygare * n1 * w3 * m3 * t2) / male.fit.bar
  
  ##### XGARe.prime #####
  XGARe.prime <- (XGARs*XGARe * n2 *w4*m1*t1 + 
                    XGARs*XGAre*.5*n2*w4*m1*t1 + 
                    XGARs*XGaRe *.5*n2*w5*m1*t1 +
                    XGARs*XGare * .5 * (r1.c*r2.c*r3.c +r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c)*n2*w5*m1*t1 + 
                    XGARs*XgARe * .5*n2*w4*m2*t1+
                    XGARs*XgAre*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3+r1*r2.c*r3.c+r1*r2*r3)*n2*w4*m2*t1 +
                    XGARs*XgaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3+r1*r2.c*r3.c + r1*r2.c*r3) * n2*w5*m2*t1 +
                    XGARs*Xgare *.5*(r1.c*r2.c*r3.c+r1*r2.c*r3.c)*n2*w5*m2*t1 +
                    XGARs*YGARe * .5 *n1.c*w4*m1*t2 + 
                    XGARs*YGAre*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3+r1*r2*r3.c + r1*r2.c*r3)*n1.c *w4*m1*t2 +
                    XGARs*YGaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3+r1*r2*r3.c+r1*r2*r3)*n1.c*w5*m1*t2+
                    XGARs*YGare*.5*(r1.c*r2.c*r3.c+r1*r2*r3.c)*n1.c*w5*m1*t2 +
                    XGARs*YgARe*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3.c + r1.c*r2*r3+r1.c*r2.c*r3)*n1.c*w4*m2*t2 +
                    XGARs*YgAre*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3)*n1.c*w4*m2*t2 +
                    XGARs*YgaRe*.5*(r1.c*r2.c*r3.c+r1.c*r2.c*r3)*n1.c*w5*m2*t2 +
                    XGARs*Ygare*.5*r1.c*r2.c*r3.c*n1.c*w5*m2*t2 +
                    # row 2
                    XGArs*XGARe*.5*n2*w4*m1*t1 +
                    XGArs*XGaRe*.5*(r1.c*r2.c*r3+r1*r2*r3+ r1*r2.c*r3+r1.c*r2*r3)*n2*w5*m1*t1+
                    XGArs*XgARe*.5*(r1.c*r2.c*r3+r1.c*r2*r3.c + r1*r2.c*r3+r1*r2*r3.c)*n2*w4*m2*t1 +
                    XGArs*XgaRe*.5*(r1.c*r2.c*r3+r1*r2.c*r3)*n2*w5*m2*t1 +
                    XGArs*YGARe*.5*(r1.c*r2.c*r3+r1*r2.c*r3.c+r1.c*r2*r3.c+r1*r2*r3)*n1.c*w4*m1*t2 +
                    XGArs*YGaRe*.5*(r1.c*r2.c*r3+r1*r2*r3)*n1.c*w5*m1*t2 +
                    XGArs*YgARe*.5*(r1.c*r2.c*r3+r1.c*r2*r3.c)*n1.c*w4*m2*t2 +
                    XGArs*YgaRe*.5*(r1.c*r2.c*r3)*n1.c*w5*m2*t2+
                    # row 3
                    XGaRs*XGARe*.5*n2*w5*m1*t1 +
                    XGaRs*XGAre*.5*(r1.c*r2*r3+r1*r2.c*r3+r1.c*r2.c*r3+r1*r2*r3)*n2*w5*m1*t1 +
                    XGaRs*XgARe*.5*(r1.c*r2*r3.c+r1.c*r2*r3+ r1*r2*r3.c+r1*r2*r3)*n2*w5*m2*t1 +
                    XGaRs*XgAre*.5*(r1.c*r2*r3+r1*r2*r3)*n2*w5*m2*t1+
                    XGaRs*YGARe*.5*(r1*r2.c*r3+r1*r2.c*r3.c+r1.c*r2*r3.c + r1.c*r2*r3)*n1.c*w5*m1*t2 +
                    XGaRs*YGAre*.5*(r1.c*r2*r3+r1*r2.c*r3)*n1.c*w5*m1*t2 +
                    XGaRs*YgARe*.5*(r1.c*r2*r3+r1.c*r2*r3.c)*n1.c*w5*m2*t2 +
                    XGaRs*YgAre*.5*r1.c*r2*r3*n1.c*w5*m2*t2 +
                    # row 4
                    XGars*XGARe*.5*(r1.c*r2.c*r3.c+r1.c*r2*r3.c+ r1*r2.c*r3.c+r1*r2*r3.c)*n2*w5*m1*t1 +
                    XGars*XgARe*.5*(r1.c*r2*r3.c+r1*r2*r3.c)*n2*w5*m2*t1 +
                    XGars*YGARe*.5*(r1*r2.c*r3.c+r1.c*r2*r3.c)*n1.c*w5*m1*t2 +
                    XGars*YgARe*.5*r1.c*r2*r3.c*n1.c*w5*m2*t2 +
                    # row 5
                    XgARs*XGARe*.5*n2*w4*m2*t1 +
                    XgARs*XGAre*.5*(r1.c*r2*r3.c+r1.c*r2.c*r3+r1*r2.c*r3+r1*r2*r3.c)*n2*w4*m2*t1 +
                    XgARs*XGaRe*.5*(r1.c*r2*r3+r1.c*r2*r3.c+r1*r2*r3+r1*r2*r3.c)*n2*w5*m2*t1+
                    XgARs*XGare*.5*(r1.c*r2*r3.c+r1*r2*r3.c)*n2*w5*m2*t1 +
                    XgARs*YGARe*.5*(r1*r2*r3.c+r1*r2.c*r3.c+r1*r2*r3+r1*r2.c*r3)*n1.c*w4*m2*t2 +
                    XgARs*YGAre*.5*(r1*r2*r3.c+r1*r2.c*r3)*n1.c*w4*m2*t2 +
                    XgARs*YGaRe*.5*(r1*r2*r3+r1*r2*r3.c)*n1.c*w5*m2*t2 +
                    XgARs*YGare*.5*r1*r2*r3.c*n1.c*w5*m2*t2 +
                    # row 6
                    XgArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XgArs*XGaRe * .5 * (r1*r2*r3 +r1.c*r2*r3) * n2 * w5 * m2 * t1 +
                    XgArs*YGARe * .5 * (r1*r2*r3 + r1*r2.c*r3.c) *n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGAre * .5 * (r1*r2.c*r3 + r1.c*r2.c*r3) * n2 * w5 *m2 *t1 +
                    XgaRs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c *w5 *m2 *t2 + 
                    XgaRs*YGAre * .5 * r1*r2.c*r3 * n1.c *w5 * m2 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1+
                    Xgars*YGARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 *t2 + 
                    # row 9
                    YGARs*XGARe * .5 * n1.c * w4 * m1 * t2 +
                    YGARs*XGAre * .5 * (r1*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 + 
                    YGARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m1 * t2 + 
                    YGARs*XGare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 + 
                    YGARs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3 +r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGARs*XgAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m1 * t2 +
                    YGArs*XGaRe * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XgARe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * (r1*r2.c*r3) * n1.c * w5 * m2 * t2 + 
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XgARe * .5 * (r1*r2*r3 + r1*r2*r3.c) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGars*XgARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 13
                    YgARs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3 + r1.c*r2*r3.c) * n1.c *w4 * m2 * t2 +
                    YgARs*XGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 + 
                    YgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 + 
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2) / female.fit.bar
  
  ##### XGAre.prime #####
  XGAre.prime <- (XGARs*XGAre * .5 * n2 * w4 * m1 * t1 +
                    XGARs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2*r3 + r1*r2.c*r3 + r1.c*r2*r3) * n2 * w5 * m1 * t1 +
                    XGARs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n2 * w4 * m2 * t1 +
                    XGARs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    XGARs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 + 
                    XGARs*Ygare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XGARe * .5 * n2 * w4 * m1 * t1 +
                    XGArs*XGAre * n2 * w4 * m1 * t1 +
                    XGArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c) * n2 * w5 * m1 * t1 +
                    XGArs*XGare * .5 * n2 * w5 * m1 * t1 +
                    XGArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XGArs*XgAre * .5 * n2 * w4 * m2 * t1 +
                    XGArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    XGArs*YGAre * .5 * n1.c * w4 * m1 * t2 +
                    XGArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1.c*r2*r3.c + r1.c*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c) * n2 * w5 * m1 * t1 +
                    XGaRs*XgAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGaRs*YGAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    XGaRs*YgAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*XGARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3 + r1.c*r2.c*r3 + r1*r2*r3) * n2 * w5 * m1 * t1 +
                    XGars*XGAre * .5 * n2 * w5 * m1 * t1 +
                    XGars*XgARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1+
                    XGars*XgAre * .5 * (r1.c*r2*r3 + r1.c*r2*r3.c + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGars*YGARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    XGars*YGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YgARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XgARs*XGare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgARs*YGAre * .5 * (r1*r2.c*r3.c +r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n2 * w4 * m2 * t1 +
                    XgArs*XGAre * .5 * n2 * w4 * m2 * t1 +
                    XgArs*XGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XgArs*XGare * .5 * (r1.c*r2*r3 + r1.c*r2*r3.c + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgArs*YGARe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 + 
                    XgArs*YGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgaRs*YGAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 + 
                    Xgars*YGARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 9
                    YGARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    YGARs*XGare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    YGARs*XgAre * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 + 
                    YGARs*Xgare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    YGArs*XGAre * .5 * n1.c * w4 * m1 * t2 +
                    YGArs*XGaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    YGArs*XGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XgARe * .5 * (r1*r2.c*r3.c +r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGaRs*XgAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XgARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 13
                    YgARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 + 
                    YgArs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1.c*r2*r3.c + r1.c*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2) / female.fit.bar
  
  ##### XGaRe.prime #####
  XGaRe.prime <- (XGARs*XGaRe * .5 * n2 * w5 * m1 * t1 +
                    XGARs*XGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m1 * t1 +
                    XGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGARs*Xgare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YGare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    XGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 + 
                    XGARs*Ygare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 + 
                    # row 2
                    XGArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m1 * t1 +
                    XGArs*XgaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGArs*YGaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    XGArs*YgaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XGARe * .5 * n2 * w5 * m1 * t1 +
                    XGaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m1 * t1 +
                    XGaRs*XGaRe * n2 * w6 * m1 * t1 +
                    XGaRs*XGare * .5 * n2 * w6 * m1 * t1 +
                    XGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGaRs*XgaRe * .5 * n2 * w6 * m2 * t1 +
                    XGaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XGaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGaRe * .5 * n1.c * w6 * m1 * t2 +
                    XGaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m1 * t2 +
                    XGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XGaRs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m1 * t1 +
                    XGars*XGaRe * .5 * n2 * w6 * m1 * t1 +
                    XGars*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGars*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n2 * w6 * m2 * t1 +
                    XGars*YGARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGars*YgARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 5
                    XgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgARs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgArs*YGaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XgaRs*XGaRe * .5 * n2 * w6 * m2 * t1 +
                    XgaRs*XGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n2 * w6 * m2 * t1   +
                    XgaRs*YGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    Xgars*YGARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 9
                    YGARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGArs*XgaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGaRe * .5 * n1.c * w6 * m1 * t2 +
                    YGaRs*XGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    YGars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m1 * t2 +
                    YGars*XgARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 + 
                    YgaRs*XGAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 + 
                    Ygars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2) / female.fit.bar
  
  ##### XGare.prime #####
  XGare.prime <- (XGARs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + 
                                        r1*r2*r3.c) * n2 * w5 * m1 * t1 + 
                    XGARs*Xgare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGARs*YGare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    XGARs*Ygare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m1 * t1 +
                    XGArs*XGare * .5 * n2 * w5 * m1 * t1 +
                    XGArs*XgaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGArs*Xgare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGArs*YGaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGArs*YGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    XGArs*YgaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m1 * t1 +
                    XGaRs*XGare * .5 * n2 * w6 * m1 * t1 +
                    XGaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*Xgare * .5 * (r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2.c*r3) * n2 * w6 * m2 * t1 +
                    XGaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGaRs*YgAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c) * n2 * w5 * m1 * t1 + 
                    XGars*XGAre * .5 * n2 * w5 * m1 * t1 +
                    XGars*XGaRe * .5 * n2 * w6 * m1 * t1 +
                    XGars*XGare * n2 * w6 * m1 * t1 +
                    XGars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XGars*Xgare * .5 * n2 * w6 * m2 * t1 +
                    XGars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 + 
                    XGars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 + 
                    XGars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m1 * t2 +
                    XGars*YGare * .5 * n1.c * w6 * m1 * t2 +
                    XGars*YgARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 + 
                    XGars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XGars*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 5 
                    XgARs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgARs*YGare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 + 
                    # row 6
                    XgArs*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgArs*YGaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 7 
                    XgaRs*XGAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XgaRs*YGAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    Xgars*XGAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2.c*r3) * n2 * w6 * m2 * t1 +
                    Xgars*XGare * .5 * n2 * w6 * m2 * t1 +
                    Xgars*YGARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 9
                    YGARs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 + 
                    YGARs*Xgare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 10 
                    YGArs*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 + 
                    YGArs*XgaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XGAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m1 * t2 +
                    YGaRs*XgAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    YGars*XGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 + 
                    YGars*XGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGars*XGare * .5 * n1.c * w6 * m1 * t2 +
                    YGars*XgARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 + 
                    YGars*XgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGars*Xgare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 13
                    YgARs*XGare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 + 
                    # row 14
                    YgArs*XGaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    Ygars*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2) / 
    female.fit.bar
  
  ##### XgARe.prime #####
  XgARe.prime <- (XGARs*XgARe * .5 * n2 * w4 * m2 * t1 + 
                    XGARs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c + r1*r2.c*r3) * n2 * w4 * m2 * t1 +
                    XGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGARs*Xgare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGARs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*YgAre * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XGARs*YgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 + 
                    XGARs*Ygare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XGArs*XgaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGArs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGars*YgARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*XGARe * .5 * n2 * w4 * m2 * t1 + 
                    XgARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgARs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgARs*XgARe * n2 * w4 * m3 * t1 +
                    XgARs*XgAre * .5 * n2 * w4 * m3 * t1 +
                    XgARs*XgaRe * .5 * n2 * w5 * m3 * t1 +
                    XgARs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgARs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 +  r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgARs*YgARe * .5 * n1.c * w4 * m3 * t2 +
                    XgARs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m3 * t2 +
                    XgARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*XGARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w4 * m2 * t1 +
                    XgArs*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3)* n2 * w5 * m2 * t1 +
                    XgArs*XgARe * .5 * n2 * w4 * m3 * t1 +
                    XgArs*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgArs*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgArs*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XgARe * .5 * n2 * w5 * m3 * t1 +
                    XgaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    Xgars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    Xgars*YGARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    Xgars*YgARe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 9 
                    YGARs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 +  r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    # row 10
                    YGArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 13
                    YgARs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgARs*XgARe * .5 * n1.c * w4 * m3 * t2 +
                    YgARs*XgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    YgARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m3 * t2 +
                    YgArs*XgaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 + 
                    YgaRs*XGAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2) / female.fit.bar
  
  ##### XgAre.prime #####
  XgAre.prime <- (XGARs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XGARs*Xgare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 *  m2 * t1 +
                    XGARs*YgAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*Ygare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w4 * m2 * t1 +
                    XGArs*XgAre * .5 * n2 * w4 * m2 * t1 +
                    XGArs*XgaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGArs*Xgare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGArs*YgARe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XGArs*YgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGaRs*YgAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGars*YgARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*XGAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w4 * m2 * t1 +
                    XgARs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgARs*XgAre * .5 * n2 * w4 * m3 * t1 +
                    XgARs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgARs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w4 * m2 * t1 +
                    XgArs*XGAre * .5 * n2 * w4 * m2 * t1 +
                    XgArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgArs*XgARe * .5 * n2 * w4 * m3 * t1 +
                    XgArs*XgAre * n2 * w4 * m3 * t1 +
                    XgArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgArs*Xgare * .5 * n2 * w5 * m3 * t1 +
                    XgArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    XgArs*YgAre * .5 * n1.c * w4 * m3 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    XgArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*XGAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XgaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgaRs*YGAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YgAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 *  m2 * t1 +
                    Xgars*XGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    Xgars*XgAre * .5 * n2 * w5 * m3 * t1 +
                    Xgars*YGARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 9
                    YGARs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGArs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XgAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 13
                    YgARs*XGAre * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 14
                    YgArs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgArs*XgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    YgArs*XgAre * .5 * n1.c * w4 * m3 * t2 +
                    YgArs*XgaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    YgArs*Xgare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2) / female.fit.bar
  
  ##### XgaRe.prime #####
  XgaRe.prime <- (XGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + 
                                        r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGARs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGARs*Ygare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGArs*YgaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*XgAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGaRs*XgaRe * .5 * n2 * w6 * m2 * t1 +
                    XGaRs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w6 * m2 * t1 +
                    XGaRs*YgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XGaRs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*XgARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XGars*YgARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 5
                    XgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgARs*XGare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgARs*XgaRe * .5 * n2 * w5 * m3 * t1 +
                    XgARs*Xgare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*XGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XgArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgArs*YGaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YgaRe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XgaRs*XGaRe * .5 * n2 * w6 * m2 * t1 +
                    XgaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XgaRs*XgARe * .5 * n2 * w5 * m3 * t1 +
                    XgaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgaRs*XgaRe * n2 * w6 * m3 * t1 + 
                    XgaRs*Xgare * .5 * n2 * w6 * m3 * t1 + 
                    XgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgaRe * .5 * n1.c * w6 * m3 * t2 +
                    XgaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m3 * t2 +
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w6 * m2 * t1 +
                    Xgars*XgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    Xgars*XgaRe * .5 * n2 * w6 * m3 * t1 + 
                    Xgars*YGARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m3 * t2 +
                    # row 9
                    YGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XgaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgAre * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgaRe * .5 * n1.c * w6 * m3 * t2 +
                    YgaRs*Xgare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m3 * t2) / 
    female.fit.bar
  
  ##### Xgare.prime #####
  Xgare.prime <- (XGARs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    XGARs*Ygare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*XgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XGArs*YgaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*XgAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    XGaRs*YgAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*XgARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XGars*XgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XGars*XgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w6 * m2 * t1 +
                    XGars*Xgare * .5 * n2 * w6 * m2 * t1 +
                    XGars*YgARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    XGars*YgAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 + 
                    XGars*YgaRe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 + 
                    XGars*Ygare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    # row 5
                    XgARs*XGare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m2 * t1 +
                    XgARs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    XgARs*YGare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    XgARs*Ygare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 + 
                    # row 6
                    XgArs*XGaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgArs*XGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n2 * w5 * m2 * t1 +
                    XgArs*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgArs*Xgare * .5 * n2 * w5 * m3 * t1 +
                    XgArs*YGaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 + 
                    XgArs*YGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 + 
                    XgArs*YgaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 + 
                    XgArs*Ygare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 + 
                    # row 7
                    XgaRs*XGAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    XgaRs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w6 * m2 * t1 +
                    XgaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n2 * w5 * m3 * t1 +
                    XgaRs*Xgare * .5 * n2 * w6 * m3 * t1 +
                    XgaRs*YGAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 + 
                    XgaRs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 + 
                    XgaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 + 
                    XgaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m3 * t2 + 
                    # row 8
                    Xgars*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n2 * w5 * m2 * t1 +
                    Xgars*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n2 * w5 * m2 * t1 +
                    Xgars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n2 * w6 * m2 * t1 +
                    Xgars*XGare * .5 * n2 * w6 * m2 * t1 +
                    Xgars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n2 * w5 * m3 * t1 +
                    Xgars*XgAre * .5 * n2 * w5 * m3 * t1 +
                    Xgars*XgaRe * .5 * n2 * w6 * m3 * t1 +
                    Xgars*Xgare * n2 * w6 * m3 * t1 +
                    Xgars*YGARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    Xgars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m3 * t2 +
                    Xgars*Ygare * .5 * n1.c * w6 * m3 * t2 +
                    # row 9
                    YGARs*Xgare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 10
                    YGArs*XgaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 + 
                    YGArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 11
                    YGaRs*XgAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 + 
                    YGaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 12
                    YGars*XgARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    YGars*XgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 + 
                    YGars*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 + 
                    YGars*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 13
                    YgARs*XGare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 + 
                    YgARs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 14
                    YgArs*XGaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 + 
                    YgArs*XgaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 + 
                    YgArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 15
                    YgaRs*XGAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 + 
                    YgaRs*XgAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 + 
                    YgaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w6 * m3 * t2 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    Ygars*XGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    Ygars*XgARe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 + 
                    Ygars*XgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m3 * t2 + 
                    Ygars*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m3 * t2 + 
                    Ygars*Xgare * .5 * n1.c * w6 * m3 * t2) / female.fit.bar 
  
  ##### YGARe.prime #####
  YGARe.prime <- (XGARs*YGARe * .5 * n1.c * w4 * m1 * t2 +
                    XGARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    XGARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YGare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGARs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*YgAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGARs*Ygare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    XGArs*YGaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGArs*YgARe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGars*YgARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 5 
                    XgARs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 8 
                    Xgars*YGARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 9
                    YGARs*XGARe * .5 * n1.c * w4 * m1 * t2 +
                    YGARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    YGARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGARs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGARs*YGARe * n1.c * w4 * m1 * t3 +
                    YGARs*YGAre * .5 * n1.c * w4 * m1 * t3 +
                    YGARs*YGaRe * .5 * n1.c * w5 * m1 * t3 +
                    YGARs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGARs*YgARe * .5 * n1.c * w4 * m2 * t3 +
                    YGARs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    YGArs*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGArs*YGARe * .5 * n1.c * w4 * m1 * t3 +
                    YGArs*YGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGArs*YgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YGArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGaRs*YGARe * .5 * n1.c * w5 * m1 * t3 +
                    YGaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 12
                    YGars*XGARe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGars*XgARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGars*YgARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 13
                    YgARs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGAre * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgARs*YGARe * .5 * n1.c * w4 * m2 * t3 +
                    YgARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 14
                    YgArs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YgArs*YGaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3) / female.fit.bar
  
  ##### YGAre.prime #####
  YGAre.prime <- (XGARs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    XGARs*YGare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YgAre * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*Ygare * .5 * (r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    XGArs*YGAre * .5 * n1.c * w4 * m1 * t2 +
                    XGArs*YGaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    XGArs*YGare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGArs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGaRs*YgAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*YGARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YgARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgArs*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*YGAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 9
                    YGARs*XGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m1 * t2 +
                    YGARs*XGare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGARs*YGAre * .5 * n1.c * w4 * m1 * t3 +
                    YGARs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGARs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w4 * m1 * t2 +
                    YGArs*XGAre * .5 * n1.c * w4 * m1 * t2 +
                    YGArs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGArs*YGARe * .5 * n1.c * w4 * m1 * t3 +
                    YGArs*YGAre * n1.c * w4 * m1 * t3 +
                    YGArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGArs*YGare * .5 * n1.c * w5 * m1 * t3 +
                    YGArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YGArs*YgAre * .5 * n1.c * w4 * m2 * t3 +
                    YGArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YGArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XGAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XGAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XgARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGars*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGars*YGAre * .5 * n1.c * w5 * m1 * t3 +
                    YGars*YgARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 13
                    YgARs*XGAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 14
                    YgArs*XGARe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgArs*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YgArs*YGAre * .5 * n1.c * w4 * m2 * t3 +
                    YgArs*YGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YgArs*YGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 15
                    YgaRs*XGAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * (r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3) / 
    female.fit.bar 
  
  ##### YGaRe.prime #####
  YGaRe.prime <- (XGARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGARs*YgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGARs*Ygare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGArs*YgaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGAre * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGaRe * .5 * n1.c * w6 * m1 * t2 +
                    XGaRs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XGaRs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*YGARe * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGars*YgARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 5
                    XgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*YGaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 8 
                    Xgars*YGARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 9
                    YGARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XGare * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGARs*XgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGARs*YGaRe * .5 * n1.c * w5 * m1 * t3 +
                    YGARs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XGaRe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t3 +
                    YGArs*YgaRe * .5 * (r1*r2*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGaRe * .5 * n1.c * w6 * m1 * t2 +
                    YGaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*YGARe * .5 * n1.c * w5 * m1 * t3 +
                    YGaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m1 * t3 +
                    YGaRs*YGaRe * n1.c * w6 * m1 * t3 +
                    YGaRs*YGare * .5 * n1.c * w6 * m1 * t3 +
                    YGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YGaRs*YgaRe * .5 * n1.c * w6 * m2 * t3 +
                    YGaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3)* n1.c * w6 * m2 * t3 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGars*XgARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    YGars*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGars*YGaRe * .5 * n1.c * w6 * m1 * t3 +
                    YGars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 14
                    YgArs*XGaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    YgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGAre * .5 * (r1*r2*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGaRe * .5 * n1.c * w6 * m2 * t3 +
                    YgaRs*YGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w6 * m2 * t3) / female.fit.bar
  
  ##### YGare.prime #####
  YGare.prime <- (XGARs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGARs*Ygare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YGaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 + # 1-r2 -> r2
                    XGArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    XGArs*YgaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YGAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGaRs*YgAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*YGARe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    XGars*YGAre * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    XGars*YGaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    XGars*YGare * .5 * n1.c * w6 * m1 * t2 +
                    XGars*YgARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    XGars*Ygare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    # row 5
                    XgARs*YGare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 6
                    XgArs*YGaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 7
                    XgaRs*YGAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    Xgars*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 9
                    YGARs*XGare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGARs*Xgare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGARs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGARs*Ygare * .5 * (r1*r2*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XGaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XGare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m1 * t2 +
                    YGArs*XgaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGArs*YGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGArs*YGare * .5 * n1.c * w5 * m1 * t3 +
                    YGArs*YgaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGArs*Ygare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XGAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGaRs*XGare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGaRs*XgAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    YGaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m1 * t3 +
                    YGaRs*YGare * .5 * n1.c * w6 * m1 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    # row 12
                    YGars*XGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m1 * t2 +
                    YGars*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m1 * t2 +
                    YGars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w6 * m1 * t2 +
                    YGars*XGare * .5 * n1.c * w6 * m1 * t2 +
                    YGars*XgARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGars*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m1 * t3 +
                    YGars*YGAre * .5 * n1.c * w5 * m1 * t3 +
                    YGars*YGaRe * .5 * n1.c * w6 * m1 * t3 +
                    YGars*YGare * n1.c * w6 * m1 * t3 +
                    YGars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YGars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    YGars*Ygare * .5 * n1.c * w6 * m2 * t3 +
                    # row 13
                    YgARs*XGare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgARs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 14
                    YgArs*XGaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgArs*YGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 15
                    YgaRs*XGAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    YgaRs*YGAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    Ygars*XGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 + 
                    Ygars*YGARe * .5 * (r1*r2*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    Ygars*YGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    Ygars*YGare * .5 * n1.c * w6 * m2 * t3) / female.fit.bar
  
  ##### YgARe.prime #####
  YgARe.prime <- (XGARs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * 
                    n1.c * w4 * m2 * t2 +
                    XGARs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGARs*Ygare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*YgARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgARs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgARs*YgARe * .5 * n1.c * w4 * m3 * t2 +
                    XgARs*YgAre * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    XgARs*YgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*YGARe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 9
                    YGARs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*XgAre * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGARs*XgaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGARs*YgARe * .5 * n1.c * w4 * m2 * t3 +
                    YGARs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YGARs*YgaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGArs*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YGArs*YgaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGaRs*YgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 12
                    YGars*XgARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 13
                    YgARs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgARs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgARs*XgARe * .5 * n1.c * w4 * m3 * t2 +
                    YgARs*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    YgARs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgARs*YGARe * .5 * n1.c * w4 * m2 * t3 +
                    YgARs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YgARs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YgARs*YgARe * n1.c * w4 * m3 * t3 +
                    YgARs*YgAre * .5 * n1.c * w4 * m3 * t3 +
                    YgARs*YgaRe * .5 * n1.c * w5 * m3 * t3 +
                    YgARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    # row 14
                    YgArs*XGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XgARe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    YgArs*XgaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgArs*YGARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YgArs*YGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgArs*YgARe * .5 * n1.c * w4 * m3 * t3 +
                    YgArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*YGARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YgARe * .5 * n1.c * w5 * m3 * t3 +
                    YgaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    Ygars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t3) / 
    female.fit.bar
  
  ##### YgAre.prime #####
  YgAre.prime <- (XGARs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    XGARs*Ygare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XGArs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XGArs*YgaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YgAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 4
                    XGars*YgARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 5
                    XgARs*YGAre * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    XgARs*YGare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgARs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    XgArs*YGaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgArs*YgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    XgArs*YgAre * .5 * n1.c * w4 * m3 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    XgArs*Ygare * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*YGAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 9
                    YGARs*XgAre * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGARs*Xgare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGARs*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XgARe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YGArs*XgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t2 +
                    YGArs*XgaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGArs*YgARe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YGArs*YgAre * .5 * n1.c * w4 * m2 * t3 +
                    YGArs*YgaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YGArs*Ygare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XgAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 12
                    YGars*XgARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 13
                    YgARs*XGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*XgAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*YGAre * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgARs*YgAre * .5 * n1.c * w4 * m3 * t3 +
                    YgARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    # row 14
                    YgArs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w4 * m2 * t2 +
                    YgArs*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w4 * m2 * t2 +
                    YgArs*XGaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgArs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w4 * m3 * t2 +
                    YgArs*XgAre * .5 * n1.c * w4 * m3 * t2 +
                    YgArs*XgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgArs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgArs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w4 * m2 * t3 +
                    YgArs*YGAre * .5 * n1.c * w4 * m2 * t3 +
                    YgArs*YGaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YgArs*YGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgArs*YgARe * .5 * n1.c * w4 * m3 * t3 +
                    YgArs*YgAre * n1.c * w4 * m3 * t3 +
                    YgArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c +  r1*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    YgArs*Ygare * .5 * n1.c * w5 * m3 * t3 +
                    # row 15
                    YgaRs*XGAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m3 * t2 +
                    YgaRs*YGAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YgaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c +  r1*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*XgAre * .5 * (r1.c*r2*r3.c + r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    Ygars*YgAre * .5 * n1.c * w5 * m3 * t3) / female.fit.bar
  
  ##### YgaRe.prime #####
  YgaRe.prime <- (XGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XGARs*Ygare * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YgaRe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGaRs*YgAre * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XGaRs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*YgARe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 5
                    XgARs*YGaRe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgARs*YGare * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*YGaRe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 7
                    XgaRs*YGARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    XgaRs*YGAre * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    XgaRs*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    XgaRs*YgARe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgAre * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    XgaRs*YgaRe * .5 * n1.c * w6 * m3 * t2 +
                    XgaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w6 * m3 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m3 * t2 +
                    # row 9
                    YGARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGARs*Xgare * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGARs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XgaRe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XgARe * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGaRs*XgAre * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGaRs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*YgARe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3)  * n1.c * w5 * m2 * t3 +
                    YGaRs*YgAre * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YGaRs*YgaRe * .5 * n1.c * w6 * m2 * t3 +
                    YGaRs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    # row 12
                    YGars*XgARe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGars*YgARe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    # row 13
                    YgARs*XGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgARs*XGare * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgARs*XgaRe * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*Xgare * .5 * (r1.c*r2*r3 + r1*r2.c*r3) * n1.c * w5 * m3 * t2 +
                    YgARs*YGaRe * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3)  * n1.c * w5 * m2 * t3 +
                    YgARs*YGare * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgARs*YgaRe * .5 * n1.c * w5 * m3 * t3 +
                    YgARs*Ygare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    # row 14
                    YgArs*XGaRe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgArs*XgaRe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgArs*YGaRe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YgArs*YgaRe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    # row 15
                    YgaRs*XGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YgaRs*XGAre * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YgaRs*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    YgaRs*XgARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgaRs*XgaRe * .5 * n1.c * w6 * m3 * t2 +
                    YgaRs*Xgare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2.c*r3) * n1.c * w6 * m3 * t2 +
                    YgaRs*YGARe * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGaRe * .5 * n1.c * w6 * m2 * t3 +
                    YgaRs*YGare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3 + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    YgaRs*YgARe * .5 * n1.c * w5 * m3 * t3 +
                    YgaRs*YgAre * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    YgaRs*YgaRe * n1.c * w6 * m3 * t3 +
                    YgaRs*Ygare * .5 * n1.c * w6 * m3 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*XgaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3 + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w6 * m3 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    Ygars*YgARe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    Ygars*YgaRe * .5 * n1.c * w6 * m3 * t3) / female.fit.bar
  
  ##### Ygare.prime #####
  Ygare.prime <- (XGARs*Ygare * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    # row 2
                    XGArs*YgaRe * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XGArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    # row 3
                    XGaRs*YgAre * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    XGaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 4
                    XGars*YgARe * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XGars*YgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    XGars*YgaRe * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    XGars*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    # row 5
                    XgARs*YGare * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    XgARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    # row 6
                    XgArs*YGaRe * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    XgArs*YGare * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    XgArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    # row 7 
                    XgaRs*YGAre * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    XgaRs*YGare * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    XgaRs*YgAre * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    XgaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m3 * t2 +
                    # row 8
                    Xgars*YGARe * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    Xgars*YGAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Xgars*YGaRe * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YGare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    Xgars*YgARe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    Xgars*YgAre * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    Xgars*YgaRe * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w6 * m3 * t2 +
                    Xgars*Ygare * .5 * n1.c * w6 * m3 * t2 +
                    # row 9
                    YGARs*Xgare * .5 * r1*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    YGARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    # row 10
                    YGArs*XgaRe * .5 * r1*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YGArs*Xgare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    YGArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YGArs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    # row 11
                    YGaRs*XgAre * .5 * r1*r2*r3 * n1.c * w5 * m2 * t2 +
                    YGaRs*Xgare * .5 * (r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGaRs*YgAre * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGaRs*Ygare * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    # row 12
                    YGars*XgARe * .5 * r1*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YGars*XgAre * .5 * (r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t2 +
                    YGars*XgaRe * .5 * (r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    YGars*Xgare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t2 +
                    YGars*YgARe * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YGars*YgAre * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YGars*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    YGars*Ygare * .5 * n1.c * w6 * m2 * t3 +
                    # row 13
                    YgARs*XGare * .5 * r1.c*r2*r3.c * n1.c * w5 * m2 * t2 +
                    YgARs*Xgare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    YgARs*YGare * .5 * (r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w5 * m2 * t3 +
                    YgARs*Ygare * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    # row 14
                    YgArs*XGaRe * .5 * r1.c*r2*r3 * n1.c * w5 * m2 * t2 +
                    YgArs*XGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m2 * t2 +
                    YgArs*XgaRe * .5 * (r1*r2.c*r3 + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgArs*Xgare * .5 * (r1*r2.c*r3.c + r1*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgArs*YGaRe * .5 * (r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgArs*YGare * .5 * (r1.c*r2*r3.c + r1.c*r2*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m2 * t3 +
                    YgArs*YgaRe * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    YgArs*Ygare * .5 * n1.c * w5 * m3 * t3 +
                    # row 15
                    YgaRs*XGAre * .5 * r1.c*r2.c*r3 * n1.c * w5 * m2 * t2 +
                    YgaRs*XGare * .5 * (r1.c*r2.c*r3 + r1.c*r2*r3.c) * n1.c * w6 * m2 * t2 +
                    YgaRs*XgAre * .5 * (r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    YgaRs*Xgare * .5 * (r1*r2.c*r3.c + r1.c*r2*r3.c + r1.c*r2.c*r3 + r1*r2*r3) * n1.c * w6 * m3 * t2 +
                    YgaRs*YGAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    YgaRs*YGare * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3.c + r1*r2*r3.c) * n1.c * w6 * m2 * t3 +
                    YgaRs*YgAre * .5 * (r1.c*r2.c*r3 + r1*r2.c*r3 + r1.c*r2*r3 + r1*r2*r3) * n1.c * w5 * m3 * t3 +
                    YgaRs*Ygare * .5 * n1.c * w6 * m3 * t3 +
                    # row 16
                    Ygars*XGARe * .5 * r1.c*r2.c*r3.c * n1.c * w5 * m2 * t2 +
                    Ygars*XGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3) * n1.c * w5 * m2 * t2 +
                    Ygars*XGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    Ygars*XGare * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1.c*r2*r3.c + r1.c*r2*r3) * n1.c * w6 * m2 * t2 +
                    Ygars*XgARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c) * n1.c * w5 * m3 * t2 +
                    Ygars*XgAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2*r3.c + r1*r2*r3) * n1.c * w5 * m3 * t2 +
                    Ygars*XgaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3 + r1*r2*r3.c) * n1.c * w6 * m3 * t2 +
                    Ygars*Xgare * .5 * n1.c * w6 * m3 * t2 +
                    Ygars*YGARe * .5 * (r1.c*r2.c*r3.c + r1*r2.c*r3.c) * n1.c * w5 * m2 * t3 +
                    Ygars*YGAre * .5 * (r1.c*r2.c*r3.c + r1.c*r2.c*r3 + r1*r2.c*r3.c + r1*r2.c*r3) * n1.c * w5 * m2 * t3 +
                    Ygars*YGaRe * .5 * (r1.c*r2.c*r3.c + r1.c*r2*r3 + r1*r2.c*r3.c + r1*r2*r3) * n1.c * w6 * m2 * t3 +
                    Ygars*YGare * .5 * n1.c * w6 * m2 * t3 +
                    Ygars*YgARe * .5 * (r1.c*r2.c*r3.c + r1*r2*r3.c + r1*r2.c*r3.c + r1.c*r2*r3.c) * n1.c * w5 * m3 * t3 +
                    Ygars*YgAre * .5 * n1.c * w5 * m3 * t3 +
                    Ygars*YgaRe * .5 * n1.c * w6 * m3 * t3 +
                    Ygars*Ygare * n1.c * w6 * m3 * t3) / female.fit.bar
  
  ##### results #####
  results <- c(XGARs.prime,
               XGArs.prime,
               XGaRs.prime,
               XGars.prime,
               XgARs.prime,
               XgArs.prime,
               XgaRs.prime,
               Xgars.prime,
               YGARs.prime,
               YGArs.prime,
               YGaRs.prime,
               YGars.prime,
               YgARs.prime,
               YgArs.prime,
               YgaRs.prime,
               Ygars.prime,
               XGARe.prime,
               XGAre.prime,
               XGaRe.prime,
               XGare.prime,
               XgARe.prime,
               XgAre.prime,
               XgaRe.prime,
               Xgare.prime,
               YGARe.prime,
               YGAre.prime,
               YGaRe.prime,
               YGare.prime,
               YgARe.prime,
               YgAre.prime,
               YgaRe.prime,
               Ygare.prime)
  
  if(return.abs.fit == T){
    abs.fit <- c(sum(results[1:16] * male.fit.bar), sum(results[17:32] * female.fit.bar))
    names(abs.fit) <- c('male' , 'female')
    stop(return(abs.fit))
  }
  
  res.names <- c('XGARs.p',
                 'XGArs.p',
                 'XGaRs.p',
                 'XGars.p',
                 'XgARs.p',
                 'XgArs.p',
                 'XgaRs.p',
                 'Xgars.p',
                 'YGARs.p',
                 'YGArs.p',
                 'YGaRs.p',
                 'YGars.p',
                 'YgARs.p',
                 'YgArs.p',
                 'YgaRs.p',
                 'Ygars.p',
                 'XGARe.p',
                 'XGAre.p',
                 'XGaRe.p',
                 'XGare.p',
                 'XgARe.p',
                 'XgAre.p',
                 'XgaRe.p',
                 'Xgare.p',
                 'YGARe.p',
                 'YGAre.p',
                 'YGaRe.p',
                 'YGare.p',
                 'YgARe.p',
                 'YgAre.p',
                 'YgaRe.p',
                 'Ygare.p')
  names(results) <- res.names
  return(results)
  
}
