#' CONSTRUCTOR for dataWeights OBJECT
#' 
#' GIVEN A \eqn{nObs x seqLen} MATRIX OF SEQUENCES AND THE INSERT OF THE SAME LENGTH
#' THIS FUNCTION COMPUTES WEIGHTS ACCORDING THE DISSIMILARITY WEIGHT MATRIX IF IT IS PROVIDED OR IN ITS ABSENCE
#' 0/1 MATCH/MISMATCH (identity matrix). 
#' IF STATS ACROSS SEQUENCES ARE TO BE GROUPED BASED ON \code{site.sets.list}, IT DOES SO ACCORDING
#' TO THE GROUPING FUNCTION SPECIFIED: \code{within.sites.acc.fn}
#' IF STATS ACROSS OBSERVATIONS ARE TO BE GROUPED BASED ON \code{obs.sequence.sets.list},
#' IT DOES SO ACCORDING TO THE GROUPING FUNCTION SPECIFIED: \code{within.obs.acc.fn}
#'
#' @param seq.mtx <to be described>
#' @param v.indices <to be described>
#' @param opt <to be described>
#'
#' @export
dataWeights <-
  function(seq.mtx, v.indices, opt) {
    
    seq.mtx [!seq.mtx %in% opt$acceptable.chars] <- NA #CHECK AGAINST THE ACCEPTABLE CHARS
    seq.length <- ncol(seq.mtx)
    obs.num <- nrow(seq.mtx)
    weights.count <- NULL
    weights <- seq.mtx
    insert.char.vector <- opt$insert.char.vector
    #~~~~~~~~~~~~SCREEN OUT SITES  BASED ON MINIMUM VARIABILITY CRITERIA~~~~~~~~~~~~~~
    if(!is.null(opt$choose.sites.func)){
      include.site <- do.call(what=opt$choose.sites.func, args=opt$choose.sites.func.arglist)
      seq.mtx <- seq.mtx[,include.site]
      seq.length <- ncol(seq.mtx)
      insert.char.vector <- insert.char.vector[include.site]
      opt$site.sets.list <- update.site.sets.list(include.site, opt$site.sets.list)
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # KL OR SKL STAT DATA OBJECT
    if(opt$stat.type == "KL" || opt$stat.type == "SKL" ){
      wData <- list(weights = weights, v.indices = v.indices, opt = opt, weights.count = weights.count) # weights.count is NULL unless SMMB option is chosen
      class(wData )<- "dataWeights"
      return(wData)
    }
    
    insert.comp.mtx <- matrix(insert.char.vector, nrow = obs.num, ncol = seq.length, byrow=T)
    rownames(insert.comp.mtx) <- rownames(seq.mtx)
    # COMPUTES WEIGHTS
    if(is.null(opt$weight.mtx)) {
      weights <- 1*(seq.mtx != insert.comp.mtx)
    } else {
      weights <- seq.mtx *NA # INITIALIZE WEIGHT
      for(pos in 1: seq.length) { #WEIGHT MTX CAN BE CONVERTED TO HASHTABLE AS KEYS could be MUCH FASTER WAY OF ACCESSING THE ELEMENTS
        expected <- insert.char.vector[pos]
        for(obs in 1: obs.num) {
          observed <- seq.mtx[obs, pos]
          weights[obs,pos] <- ifelse(is.na(observed), NA, opt$weight.mtx[expected, observed])
        }
      }   
    }
    
    
    # COMPUTES weights.count FOR mimic.smmb 
    if(opt$stat.type ==  "mimic.smmb"){
      weights.count <- 1*!is.na(weights)
    }
    
    #ACCUMULATION ACROSS OBSERVATION
    if(!is.null(opt$obs.vSeq.sets.list) || !is.null(opt$obs.pSeq.sets.list)){
      wData <- list(weights = weights, v.indices = v.indices, opt = opt, weights.count = weights.count)
      w.obs.acc <- obs.accumulate(weight.obj = wData)
      weights <- w.obs.acc$weights
      v.indices <- w.obs.acc$v.indices
      weights.count <- w.obs.acc$weights.count
    } 
    
    # ACCUMULATES WEIGHTS ACROSS SITES IN THE SITE SETS
    if(!is.null(opt$site.sets.list)) {
      obs.names <- rownames(weights)
      site.sets.names <- names(opt$site.sets.list)
      if(is.null(site.sets.names))site.sets.names <- paste("site.set",1:length(opt$site.sets.list),sep = "")
      obs.num <- nrow(weights)
      weights<-
        matrix(
          unlist(
            lapply(opt$site.sets.list, FUN=function(set) apply(matrix(weights[,set],nrow= obs.num), MARGIN=1, FUN= opt$within.sites.acc.fn))
          ),
          nrow=nrow(weights),
          byrow = F
        )
      if(opt$stat.type ==  "mimic.smmb"){
        weights.count <-
          matrix(
            unlist(
              lapply(opt$site.sets.list, FUN=function(set) apply(matrix(weights.count[,set],nrow= obs.num), MARGIN=1, FUN= opt$within.sites.acc.fn))
            ),
            nrow=nrow(weights),
            byrow = F
          )
        rownames(weights.count) <- obs.names
        colnames(weights.count) <- site.sets.names 
      }
      
      rownames(weights) <- obs.names
      colnames(weights) <- site.sets.names 
      
    }
    
    wData <- list(weights = weights, v.indices = v.indices, opt = opt, weights.count = weights.count) # weights.count is NULL unless SMMB option is chosen
    
    class(wData )<- "dataWeights"
    return(wData)
  }

#' This function updates and returns variables of an object of the class dataWeights
#' @param dw an object of class dw to be modified
#' @param ... comma separated named variables to replace the variables of the same name in dw. Variables with invalid names are ignored.
update.dataWeights <- function(dw,...){
  update.list <- list(...)
  valid.names <- names(update.list)[names(update.list) %in% names(dw)]
  for(name in valid.names){
    dw[[name]] <- update.list[[name]]
  } 
  return(dw)
}
