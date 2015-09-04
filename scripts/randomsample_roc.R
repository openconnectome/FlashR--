require(Matrix)

con <- url("http://www.cis.jhu.edu/~parky/MRN/cci.txt")
cci <- scan(con)
close(con)
m <- length(cci)
cci.na <- which(is.na(cci))

Deltamat <- matrix(0,nrow=m,ncol=m)
for(i in 1:m) for(j in 1:m) Deltamat[i,j] <- abs(cci[i]-cci[j])
diag(Deltamat) <- NA

con <- url("http://www.cis.jhu.edu/~parky/MRN/fibergraph.Rbin")
load(con)
close(con)
m <- length(fibergraph.list)
n <- nrow(fibergraph.list[[1]])

A.list <- list()
for(i in 1:m) {
  gi <- fibergraph.list[[i]] ## upper-triangle only
  gi <- (gi + t(gi))
  gi.log <- log(gi + 1)
  A.list[[i]] <- gi.log + Diagonal(x=rowSums(gi.log))/(n-1) # diagonal augmentation
}





OR=T
AND=F
Tvmat.list <- list()
for(v in 1:n)
{
  Tvmat.list[[v]] <- matrix(0,m,m)
  for(i in 1:m)
    for(j in 1:m)
    {
      if(AND)  Tvmat.list[[v]][i,j] <- norm(as.matrix(A.list[[i]][v,v] - A.list[[j]][v,v]) , type="F")
      if(OR)   Tvmat.list[[v]][i,j] <- norm(as.matrix(A.list[[i]][v, ] - A.list[[j]][v, ]) , type="F")
    }
}

vslope = NULL
for(v in 1:n)
  vslope[v] <- lm(c(Tvmat.list[[v]][upper.tri(Tvmat.list[[v]])])~c(Deltamat[upper.tri(Deltamat)]))$coeff[2]
# > order(vslope)
#  [1] 33 49  6 41 43 15 34 30 51 45 58 66 32 57 22 69  9 17 12 23 27 31 67 70 13 37 50 42 47 40
# [31] 14 26 24  7 62 36 59  1 19 63  5 65 68 61  2 53 21 39 52 46  3 44 28 10 11  8 48 38 20 18
# [61] 55 29 16 35  4 64 54 56 25 60


nmc <- 10000
vpval = NULL
for(v in 1:n)
{
  set.seed(12345)
  permslope <- NULL
  for (mc in 1:nmc) {
    sss <- sample(1:nrow(Tvmat.list[[v]]))
    Zmat <- Tvmat.list[[v]][sss,sss]
    permslope[mc] <- lm(c(Zmat[upper.tri(Zmat)])~c(Deltamat[upper.tri(Deltamat)]))$coefficients[2]
  }
  vpval[v] <- sum(permslope>vslope[v])/nmc
}
# > order(vpval)
#  [1] 60 56 25 29 54 52  4 48 64 38 46 35  8 28 11  3 18 62 16 21 61 39  1 10  5 36 55 63 53 20
# [31] 24  7  2 14 26 19 44 59 42 40 13 65 50 68 27 47 70 37 31 67 23 69 12 17  9 22 57 66 32 58
# [61] 45 51 30 43 15 34 41  6 49 33

#> cor(vpval,vslope)
#[1] -0.9192301







W = 1:70 ; OR=T ; AND=F
Tmat <- matrix(0,m,m)
for(i in 1:m)
  for(j in 1:m)
  {
    if(AND)  Tmat[i,j] <- norm(A.list[[i]][W,W] - A.list[[j]][W,W] , type="F")
    if(OR)   Tmat[i,j] <- norm(A.list[[i]][W, ] - A.list[[j]][W, ] , type="F")
  }
ccilm <- lm(c(Tmat[upper.tri(Tmat)])~c(Deltamat[upper.tri(Deltamat)]))
slopev <- ccilm$coeff[2]
nmc <- 10000
set.seed(12345)
permslope <- NULL
for (mc in 1:nmc) {
  sss <- sample(1:nrow(Tmat))
  Zmat <- Tmat[sss,sss]
  permslope[mc] <- lm(c(Zmat[upper.tri(Zmat)])~c(Deltamat[upper.tri(Deltamat)]))$coefficients[2]
}
pval <- sum(permslope>slopev)/nmc
# W =  1:70 ; OR=T ; AND=F  =>  slope = 0.2412 ; pval = 0.0108
# W =  1:35 ; OR=T ; AND=F  =>  slope = 0.1655 ; pval = 0.0168
# W = 36:70 ; OR=T ; AND=F  =>  slope = 0.1754 ; pval = 0.0111
# W =  1:70 ; OR=F ; AND=T  =>  slope = 0.2412 ; pval = 0.0108
# W =  1:35 ; OR=F ; AND=T  =>  slope = 0.0658 ; pval = 0.0687
# W = 36:70 ; OR=F ; AND=T  =>  slope = 0.0759 ; pval = 0.0416
# W = order(vpval)[ 1:35]  ; OR=F ; AND=T  =>  slope =  0.2870 ; pval = 0.0000
# W = order(vpval)[36:70]  ; OR=F ; AND=T  =>  slope = -0.0157 ; pval = 0.5935
# W = order(vslope)[ 1:35] ; OR=F ; AND=T  =>  slope = -0.0181 ; pval = 0.6220
# W = order(vslope)[36:70] ; OR=F ; AND=T  =>  slope =  0.3068 ; pval = 0.0000
# W = order(vslope)[61:70] ; OR=F ; AND=T  =>  slope =  0.1160 ; pval = 0.0002





X1 <- X2 <- matrix(0,n,n)
for(i in 1:m)
{
  X1 <- X1 + A.list[[i]]
  X2 <- X2 + A.list[[i]]^2
}
Abar <- X1/m
Avar <- X2/m - (X1/m)^2
Asd <- sqrt(Avar)

print(date())
nmc <- 10
mm <- 109 # 109 !!!!! (now) aot 114 ?????  setdiff(1:114,cci.na)
K <- 20
uuu <- 1
roc <- matrix(0,nmc,n)
for(mc in 1:nmc)
{
  set.seed(mc+12345)
  Wstar <- sample(n,K,prob=rank(vslope))
  #Wstar <- sample(n,K)
  #Wstar = order(vslope)[51:70]
  
  AA.list <- list()
  for(i in 1:m)
  {
    AA.list[[i]] <- Abar + matrix(rnorm(n^2,0,uuu*as.vector(Asd)),n,n)
    AA.list[[i]][Wstar,Wstar] <- A.list[[i]][Wstar,Wstar]
  }
  
  thesem <- setdiff(1:114,cci.na)
  thesem <- sample(thesem,mm)
  Deltamatthesem <- Deltamat[thesem,thesem]
  
  OR=T
  AND=F
  thisTvmat.list <- list()
  for(v in 1:n)
  {
    thisTvmat.list[[v]] <- matrix(0,mm,mm)
    for(thisi in 1:mm)
      for(thisj in 1:mm)
      {
        i = thesem[thisi]
        j = thesem[thisj]
        if(AND)  thisTvmat.list[[v]][thisi,thisj] <- norm(as.matrix(AA.list[[i]][v,v] - AA.list[[j]][v,v]) , type="F")
        if(OR)   thisTvmat.list[[v]][thisi,thisj] <- norm(as.matrix(AA.list[[i]][v, ] - AA.list[[j]][v, ]) , type="F")
      }
  }
  
  thisvslope <- NULL
  for(v in 1:n)
    thisvslope[v] <- lm(c(thisTvmat.list[[v]][upper.tri(thisTvmat.list[[v]])])~c(Deltamatthesem[upper.tri(Deltamatthesem)]))$coeff[2]
  
  roc[mc,] <- cumsum(order(thisvslope,decreasing=T) %in% Wstar)
}
print(date())
plot((1:70-apply(roc,2,mean))/(n-K),apply(roc,2,mean)/K,ylim=c(0,1),xlab="P[FD]",ylab="P[TD]",type="l",lwd=3,col="black")
segments(0,0,1,1,lty=3)
for(mc in 1:nmc) points((1:70-apply(roc,2,mean))/(n-K),roc[mc,]/K,type="l",lty=2,lwd=1,col="gray")
points((1:70-apply(roc,2,mean))/(n-K), apply(roc/K,2,mean) + apply(roc/K,2,sd)/sqrt(nmc) ,ylim=c(0,1),xlab="P[FD]",ylab="P[TD]",type="l",lwd=1,col="black")
points((1:70-apply(roc,2,mean))/(n-K), apply(roc/K,2,mean) - apply(roc/K,2,sd)/sqrt(nmc) ,ylim=c(0,1),xlab="P[FD]",ylab="P[TD]",type="l",lwd=1,col="black")
#> abline(h=0.5,col="red")
#> abline(v=0.1,col="red")
# Rplot.m109.pdf  (nmc=100 => 12 hours 42 minutes)  [nmc=10  => 83 minutes]
# Rplot.m50.pdf   (nmc=100 =>  3 hours           )
# Rplot.m10.pdf   (nmc=100 =>           8 minutes)