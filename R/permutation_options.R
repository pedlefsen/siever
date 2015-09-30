#' This function generates an object of class permutation.options containing a valid set of permutation options
#' @param perm.num                  Number of permutations
#' @param refined.perm.factor       <to be described>
#' @param p.value.refined.threshold <to be described>
#' @param with.replacement          <to be described>
#' @param seed                      Random seed 
#' @export
#' 
permutation.options <- function(perm.num = 0, 
                                refined.perm.factor = 1,
                                p.value.refined.threshold = .15,
                                with.replacement = F,
                                seed = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INPUT VALIDAITON~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  numeric.params <- c("perm.num","refined.perm.factor","p.value.refined.threshold")
  are.not.numeric <- c(!is.numeric(perm.num), 
                       !is.numeric(refined.perm.factor),
                       !is.numeric(p.value.refined.threshold))
  
  if(any(are.not.numeric))
    stop(paste(numeric.params[are.not.numeric],"need to be numeric"))
  
  if(!is.null(seed) & !is.numeric(seed))
    stop("seed needs to be NULL or a valid integer")
  
  if(perm.num!=0 & !is.natural(perm.num)) 
    stop("perm.num needs to be zero or a natural number")
  
  if(refined.perm.factor<1)
    stop("refined.perm.factor must be 1 or larger")
  
  if(p.value.refined.threshold <0 || 1 < p.value.refined.threshold ) 
    stop("p.value.refined.threshold needs to be between 0 and 1")
  
  if(!is.logical(with.replacement))
    stop("with.replacement needs to be logical")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  refined.perm.num <- round((refined.perm.factor - 1)) * perm.num
  
  perm.opt <- list(perm.num = perm.num, 
                       refined.perm.num = refined.perm.num ,
                       p.value.refined.threshold = p.value.refined.threshold,
                       with.replacement = with.replacement,
                       seed = seed)
  class(perm.opt) <- "permutation.options"
  return(perm.opt)
}
                            
