#' THIS FUNCTION FINDS FINDS THE SITES MEETING THE MINIMUM VARIABILITY CRITERIA AND ADJUST THE SITES LIST ACCORDINGLY. 
#'  object in the current sieveR package. 
#' 
#' @param data.mtx              <to be described>
#' @param insert.char.vector    <to be described>
#' @param obs.vSeq.sets.list    <to be described>
#' @param obs.pSeq.sets.list    <to be described>
#' @param min.mismatches.to.insert <to be described>
#' @param min.matches.to.insert    <to be described>
#' @param alwaysScreenInsertGaps   <to be described>
#' @export
get.variable.sites  <- function(data.mtx,
                                insert.char.vector,
                                obs.vSeq.sets.list,
                                obs.pSeq.sets.list,
                                min.mismatches.to.insert = 3, 
                                min.matches.to.insert = min.mismatches.to.insert, 
                                alwaysScreenInsertGaps = FALSE){
  
  all.obs.num <- nrow(data.mtx)
  seq.length <- ncol(data.mtx)
  all.obs.sets <- c(obs.vSeq.sets.list, obs.pSeq.sets.list)
  all.obs.names <- rownames(data.mtx)
  insert.comp.mtx <- matrix(insert.char.vector, nrow = all.obs.num, ncol = seq.length, byrow=T)
  match.mtx <- 1*(data.mtx == insert.comp.mtx)
  mismatch.mtx <- 1 - match.mtx
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #SETS THE NA TO 0 IN BOTH MATCH AND MISMTACH MATRIX
  match.mtx[is.na(match.mtx)] <- 0 
  mismatch.mtx[is.na(mismatch.mtx)] <- 0
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  sets.indices.mtx <- diag(x=1,nrow=all.obs.num,ncol=all.obs.num)
  if(!is.null(obs.vSeq.sets.list) || !is.null(obs.vSeq.sets.list)){
    gen.01.row <- function (len,indices){
      binary.row <- rep(0,len) 
      binary.row[indices] <- 1 
      return(binary.row)
    } 
    sets.indices.mtx.flat <- sapply(all.obs.sets, FUN = function(.set){
      obs.set.indices <- sapply(.set,FUN= function(obs) which(all.obs.names == obs),simplify=T,USE.NAMES=F)
      return(gen.01.row(len=all.obs.num,indices = obs.set.indices))
    }, simplify = T, USE.NAMES = F)
    
    sets.indices.mtx <- matrix(sets.indices.mtx.flat,ncol = all.obs.num, byrow = T)
  }
  
  
  num.matches.sites <- apply(sets.indices.mtx %*% match.mtx > 0, MARGIN=2,FUN=sum)
  num.mismatches.sites <- apply(sets.indices.mtx %*% mismatch.mtx > 0, MARGIN=2,FUN=sum)
  
  include.site <- min.matches.to.insert <= num.matches.sites & min.mismatches.to.insert <= num.mismatches.sites
  
  if(alwaysScreenInsertGaps){
    include.site[insert.char.vector == "-"] <- FALSE
  } 
  return(include.site)
}