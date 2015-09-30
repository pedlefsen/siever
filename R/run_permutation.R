#' This is an internal function which handle permutation (with or w/o refinement)
#' @param dw data weight object
#' @param stats test statistics
#' @param perm.opt permutation.options object
#' @param get.stat.func function for computing the statistics
 
run.permutation <- function (dw,stats,perm.opt,get.stat.func){
  
  sites.names <- names(stats$test.stats)
  seq.length <- ncol(dw$weights)
  vaccinee.num <- ifelse(is.null(dw$opt$obs.vSeq.sets.list), length(dw$v.indices),length(dw$opt$obs.vSeq.sets.list)) # TODO RENAME num.vaccinees
  placebo.num <-  ifelse(is.null(dw$opt$obs.pSeq.sets.list), nrow(dw$weights) - length(dw$v.indices),  length(dw$opt$obs.pSeq.sets.list)) # TODO RENAME num.placeboees
  total.seq.num <- vaccinee.num + placebo.num # TODO RENAME num.seqs
  

  
  perm.stats <- matrix(NA, nrow = max(perm.opt$perm.num,1), ncol = seq.length,dimnames=list(NULL, sites.names))
  
  if(!is.null(perm.opt$seed))
    set.seed(perm.opt$seed)
  
  for(i in 1:perm.opt$perm.num){
    vac.smp <- sample(x=1:total.seq.num,size=vaccinee.num,replace=perm.opt$with.replacement)
    perm.stats[i,] <- get.stat.func(weight.obj=dw,
                               (1: total.seq.num)[vac.smp],
                               (1 : total.seq.num)[-vac.smp],
                               perm = T)$test.stats
  }
  p.values <- get.p.value(stats = stats, perm.stats = perm.stats ,type = dw$opt$stat.type)
  
  #~~~~~~~~~~refinement~~~~~~~~~~~~~~~
  if(perm.opt$refined.perm.num > 1){
    include.site.refined <- p.values < perm.opt$p.value.refined.threshold
    refined.seq.num <- sum(include.site.refined)
    
    if(refined.seq.num > 0){
      
      dw <- update(dw,weights = dw$weights[,include.site.refined,drop = F])
      #site.names <- names(stats$test.stats[include.site.refined])
      stats <- lapply(stats,FUN=function(.stat) return(.stat[include.site.refined]))
      
      refined.perm.stats <- matrix(NA, nrow = perm.opt$refined.perm.num, ncol = refined.seq.num, dimnames=list(NULL,names(stats$test.stats)))
      for(i in 1:perm.opt$refined.perm.num){
        vac.smp <- sample(x=1:total.seq.num,size=vaccinee.num,replace=perm.opt$with.replacement)
        refined.perm.stats[i,] <- get.stat.func(weight.obj=dw,
                                                (1: total.seq.num)[vac.smp],
                                                (1 : total.seq.num)[-vac.smp],
                                                perm = T)$test.stats
      }
      perm.stats <- rbind(perm.stats[,include.site.refined,drop=F],refined.perm.stats)
      p.values <- get.p.value(stats = list(test.stats = stats$test.stats), perm.stats = perm.stats ,type = dw$opt$stat.type)
    }
    else
      warning("None of the sites passed the refinement threshold criteria. Total of ",
              perm.opt$perm.num," permutations were completed!")
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return(list(stats=stats, perm.stats = perm.stats, p.values = p.values))
}