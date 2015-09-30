#' Main GWJ sieve analysis driver
#' 
#' GWJsieve serves AS THE PRIMARY FUNCTION FOR PERFORMING VARIOUS
#' ANALYSES USING tTEST and/or fTEST WITH OR WITHOUT PERMUTATIONS
#'
#' @param data.mtx  A \eqn{(v+p) x seqlen}  MATRIX HAVING \eqn{v} SEQUENCES FROM VACCINE GROUP AND \eqn{p} FROM PLACEBO GROUP
#' @param vaccine.indices  Indices OF THE ROWS OF \code{data.mtx} CORRESPONDING TO THE VACCINEE COHORT,
#'        THE COMPLEMENT IS ASSUMED TO BE FROM THE PLACEBO COHORT
#' @param perm.opt An object of type permutation.options containing the valid permutation parameters
#' @param wOpt AN OBJECT OF TYPE \code{\link{weightopt}} WHICH CONTAINS THE OPTIONS FOR TRANSFORMING DATA IN THE WEIGHT MATRIX
#' @param sw.test this is temporary option to allow for software testing (to be modified)
#' @export
#' @keywords sieve-analysis T-statistic GWJ
#'
GWJsieve <-
  function(data.mtx, vaccine.indices, perm.opt = permutation.options(), wOpt, sw.test = F) {
    
    if(class(perm.opt)!= "permutation.options")
      stop("perm.opt needs to be a permutation.control object")
    
    if(class(wOpt)!= "weightopt")
      stop("wOpt needs to be a weightopt object")
    
    if(!sw.test & perm.opt$perm.num == 0 & wOpt$stat.type != "tStat") 
      stop(paste(wOpt$stat.type,"requires perm.num > 0"))
   
    

    get.stat.funcs <- list(tStat = GWJ_tTest.alt, mimic.smmb = smmb.test, KL = klStat.alt, SKL = klStat.alt)
    get.stat <- get.stat.funcs[[wOpt$stat.type]]

    dw <- dataWeights(seq.mtx = data.mtx, v.indices= vaccine.indices, opt = wOpt)
    
    
    global.test <- !is.null(wOpt$stats.across.site.sets.accumulation.fns)
    global.stats <- NA
    global.p.value <- NA
    
    stats <- get.stat(weight.obj = dw,
                      vaccine.indices = dw$v.indices,
                      placebo.indices = (1:nrow(dw$weights)) [-dw$v.indices],
                      perm = F)
    if(global.test){
      global.stats <- unlist(lapply(wOpt$stats.across.site.sets.accumulation.fn, function(.fn) .fn(stats$test.stats)))
    }
    if(sw.test) return(stats$test.stats)#TO BE DELETED!
    
    
    if(perm.opt$perm.num == 0){
      p.values <- get.p.value(stats = stats, type = wOpt$stat.type)
    }
    else{
      perm.res <- run.permutation(dw=dw,stats = stats, perm.opt=perm.opt,get.stat.func=get.stat)
      p.values <- perm.res$p.values
      stats <- perm.res$stats
      
      if(global.test){
        site.cum.fns <- wOpt$stats.across.site.sets.accumulation.fns
        global.perm.stats.flat <- unlist(lapply(site.cum.fns, function(.fn) apply(perm.res$perm.stats,MARGIN=1,FUN=.fn)))
        global.perm.stats <- matrix(global.perm.stats.flat,nrow=nrow(perm.res$perm.stats), byrow=F,dimnames=list(rownames(perm.res$perm.stats), names(site.cum.fns)))
        global.p.value <- get.p.value(stats = list(test.stats = global.stats), perm.stats = global.perm.stats ,type = wOpt$stat.type)
      }
    } 
    
    
    analysis.results <- list(per.site.p.value = p.values, 
                             per.site.test.stats = stats$test.stats,
                             global.p.values = global.p.value,
                             global.test.stats = global.stats,
                             opt = wOpt)
    return(analysis.results)
  }