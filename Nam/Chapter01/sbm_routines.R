#' Testing Routine for SBM 
#'
#' @param Gdata a single igraph object
#' @param maxK The upper bound on the number of components
#' @param nboot The number of bootstrap replicates to be used 
#' @return
#' 
sbm.test <- function(Gdata,maxK=3,nboot=10) {
    g.adj = as.matrix(get.adjacency(Gdata))
    nvertex = nrow(g.adj)
    foreach(nK=1:maxK,.errorhandling = 'remove') %do% {
        fit.null = do.fit(Gdata,g.adj,nvertex,nK =nK) 
        fit.alt = do.fit(Gdata,g.adj,nvertex,nK =nK+1) 
        my.stat = fit.alt$negloglik - fit.null$negloglik                                                  

        boot.dist.null = do.boot(numboot=nboot,nvertex,fit.null,nullK=nK,altK=nK+1)       
                        
        my.ecdf = ecdf(boot.dist.null)
        list(nK=nK, fitted=fit.null, stat=my.stat, pval= my.ecdf(my.stat))
    }    
}

#' Power Calc. Routine for SBM 
#'
#' @param nullparam 
#' @param altparam
#' @param alpha
#' @param maxmc
#' @return
#' 
sbm.power <- function(nvertex,nullK, altK, nullparam, altparam, alpha=0.05, maxmc = 10) {
    
    boot.dist.null = do.boot(numboot=maxmc, nvertex, nullparam, nullK=nullK, altK=altK,)
    boot.dist.alt  = do.boot(numboot=maxmc, nvertex, altparam,  nullK=nullK, altK=altK)
    
    cval = quantile(boot.dist.null,alpha)
    ecdf(boot.dist.alt)(cval)
}

do.fit <- function(g,g.adj,nvertex,nK){                                                                                   
    g.ase = adjacency.spectral.embedding(g,nK)                                                                            
    g.ase.c = Mclust(g.ase$X,nK)                                                                                          
    
    VCop = foreach(itr=1:nK,.combine='rbind') %do% {                                                                      
        retval = rep(0,nvertex)                                                                                           
        retval[which(g.ase.c$classification==itr)]=1                                                                
        retval                                                                                                      
    }                                                                                                               
    VCop = rbind(VCop)                                                                                              
    
    normalizer = cbind(rowSums(VCop))                                                                               
    normalizer = normalizer %*% t(normalizer)                                                                       
    memb = foreach(r.itr=1:nK) %do% which(VCop[r.itr,]>0)                                                           
    
    g.mean.b = VCop %*% g.adj %*% t(VCop) / normalizer                                                              
    g.mean.full = matrix(0,nvertex,nvertex)                                                                         
    for(itr.r in 1:nK) {                                                                                            
        for(itr.c in 1:nK) {                                                                                        
            g.mean.full[memb[[itr.r]],memb[[itr.c]]] = g.mean.b[itr.r,itr.c]                                        
        }                                                                                                           
    }                                                                                                               
    
    negloglik = -sum(log(g.mean.full^g.adj)[upper.tri(g.mean.full, diag = TRUE)])                                   
    
    list(g.mean.b=g.mean.b,g.mean.full=g.mean.full,negloglik=negloglik,bstruct=sapply(memb,length))                 
}

do.boot <-function(numboot=1,nvertex,model,nullK=2,altK=3, ndopar=1,myerrorhandling='remove') {                      
    if(ndopar>1) {                                                                                                  
        foreach(boot.itr = 1:numboot,.combine='c',.errorhandling = 'remove') %dopar% {                              
            g.boot.null <- sbm.game(nvertex, pref.matrix=model$g.mean.b, block.sizes=model$bstruct)                 
            g.adj.boot.null = as.matrix(get.adjacency(g.boot.null))                                                 
            fit.null.boot = do.fit(g.boot.null,g.adj.boot.null,nvertex,nK=nullK)                                    
            fit.alt.boot = do.fit(g.boot.null,g.adj.boot.null,nvertex,nK=altK)                                      
            (fit.alt.boot$negloglik - fit.null.boot$negloglik)                                                      
        }                                                                                                           
    } else {                                                                                                        
        foreach(boot.itr = 1:numboot,.combine='c',.errorhandling = 'remove') %dopar% {                                 
            
            g.boot.null <- sbm.game(nvertex, pref.matrix=model$g.mean.b, block.sizes=model$bstruct)                 
            g.adj.boot.null = as.matrix(get.adjacency(g.boot.null))                                                 
            fit.null.boot = do.fit(g.boot.null,g.adj.boot.null,nvertex,nK=nullK)                                    
            fit.alt.boot = do.fit(g.boot.null,g.adj.boot.null,nvertex,nK=altK)                                      
            (fit.alt.boot$negloglik - fit.null.boot$negloglik)                                                      
        }                                                                                                           
    }
}
  