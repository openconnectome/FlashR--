set.seed(123)


require(foreach)
require(doMC)
registerDoMC(23)

numr = 20
numc = 1
my.p = 0
retval = foreach(mcitr=1:100,.combine='rbind',.errorhandling = 'remove') %do% {
    Z = matrix(rnorm(numr*numc,sd = 0.1),numr,numc)
    Y = 1.0*(Z > 0) + matrix(sample(0:1,numc*numr,prob=c(1-my.p,my.p),replace=TRUE),numr,numc)
    Y2 = matrix(rnorm(numr*numc),numr,numc)
    
    A = myfun(Z)
    B = myfun(Y)
    B2 = myfun(Y2)
    
    retval1 = mean(A * B)/sqrt(mean(A*A) * mean(B*B))
    retval2 = mean(A * B2)/sqrt(mean(A*A) * mean(B2*B2))    
    c(retval1,retval2)
}


hist(retval[,1],col=rgb(0,0,1,0.1),xlim=c(0,1))
hist(retval[,2],col=rgb(1,0,0,0.1),xlim=c(0,1),add=TRUE)

