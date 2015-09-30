#' weightopt CLASS CONSTRUCTOR
#'
#' Note: we may want to reconsider the structure of this constructor
#' to make it more "polymorphic".  E.g. Replace the fixed arguments
#' with "..."
#'
#' @param weight.mtx <to be described>
#' @param stat.type text string describing the types of test statistics to report.  
#' @param insert.char.vector <to be described>
#' @param site.sets.list <to be described>
#' @param obs.vSeq.sets.list <to be described>
#' @param obs.pSeq.sets.list <to be described>
#' @param acceptable.chars <to be described>
#' @param within.sites.acc.fn <to be described>
#' @param within.obs.acc.fn <to be described>
#' @param choose.sites.func <to be described> 
#' @param choose.sites.func.arglist <to be described>
#' @param stats.across.site.sets.accumulation.fns <to be described>
#' @export
weightopt <-
  function (weight.mtx=NULL,
            stat.type = c("tStat","mimic.smmb","KL","SKL"),
            insert.char.vector = NULL,
            site.sets.list = NULL,
            obs.vSeq.sets.list = NULL,
            obs.pSeq.sets.list = NULL,
            acceptable.chars = NULL,
            within.sites.acc.fn = sum,
            within.obs.acc.fn = mean,
            choose.sites.func = NULL,
            choose.sites.func.arglist = NULL,
            stats.across.site.sets.accumulation.fns = NULL ){
    
    stats.req.insert <- c("tStat", "mimic.smmb")
    #~~~~~~~~~~~~~~~~~~~~~~~VALIDATING THE INPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    stat.type <- tryCatch( match.arg(stat.type),error = function(e){
      message("test.type needs to be tStat, mimic.smmb, KL or SKL")
      return(NA)
    },finally={})
    if(is.na(stat.type)){
      return()
    }
    
    if(any(stat.type == stats.req.insert) & is.null(insert.char.vector)){
      stop(paste("Error:", stat.type, "requires valid insert.char.vector"),call.=F)
    }
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mimic.smmb <- stat.type == "mimic.smmb"
    KL <- stat.type == "KL"
    SKL <- stat.type == "SKL"
    if(KL || SKL){
      site.sets.list <- NULL
      obs.vSeq.sets.list <- NULL
      obs.pSeq.sets.list <- NULL
      message("NOTE: KL and SKL option currently does not support accumulation within sites or observation sets")
    }
    
    if( is.null( site.sets.list ) ) {
      site.sets.have.cardinality.one <- TRUE
    } else {
      site.sets.have.cardinality.one <-
        all(unlist(lapply(site.sets.list, FUN=function(set)length(set)== 1)))
    } 
    
    if( is.null( weight.mtx ) ){ 
      if(is.null(acceptable.chars)) {
        acceptable.chars <- c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" )
      }
    } else { 
      acceptable.chars <- setdiff( rownames( weight.mtx ), 'X' )
      message("NOTE: When weight matrix is provided acceptable characters are overridden based on row names of the matrix\n") 
    }
    
    if( mimic.smmb ) {
      #within.obs.acc.fn <- function(nums) return(mean(nums, na.rm=T))
      within.obs.acc.fn <- sum
      message("NOTE: When mimic option is chosen function for accumulation-across-observations-in-a-group is always overridden to sum\n") 
    }
    
    opts <-
      list(
        weight.mtx = weight.mtx,
        stat.type = stat.type,
        insert.char.vector = insert.char.vector,
        site.sets.list = site.sets.list,
        site.sets.have.cardinality.one = site.sets.have.cardinality.one,
        obs.vSeq.sets.list = obs.vSeq.sets.list,
        obs.pSeq.sets.list = obs.pSeq.sets.list,
        acceptable.chars = acceptable.chars,
        within.sites.acc.fn = within.sites.acc.fn,
        within.obs.acc.fn = within.obs.acc.fn,
        choose.sites.func = choose.sites.func,
        choose.sites.func.arglist = choose.sites.func.arglist,
        stats.across.site.sets.accumulation.fns = stats.across.site.sets.accumulation.fns,
        call = match.call()
      )
    
    class(opts) <- "weightopt"
    return(opts)
  }

#' This function updates and returns variables of an object of the class weightopt
#' @param wOpt an object of class weightopt to be modified
#' @param ... comma separated named variables to replace the variables of the same name in wOpt. Variables with invalid names are ignored.
update.weightopt <- function(wOpt,...){
  update.list <- list(...)
  valid.names <- names(update.list)[names(update.list) %in% names(wOpt)]
  for(name in valid.names){
    wOpt[[name]] <- update.list[[name]]
  } 
  return(wOpt)
}