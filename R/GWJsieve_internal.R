#' Convert an environment full of sequences into a matrix full of sequences
#'
#' This function converts \code{seq.env}, an environment containing \eqn{n} observations of sequences of length \eqn{m},
#' into an \eqn{n x m} matrix. It requires the sequences to have the same length. Otherwise, it removes all
#' but sequences of the 1st mode length, in which case it generates a warning.
#'
#' @param seq.environment environment containing aligned sequence strings.
#' @return matrix, rows are sequences, columns are positions, row labels are sequence names
convert.env.to.matrix <-
    function (seq.environment){
       obs.num <- length(unlist(as.list(seq.environment)))
       seq.lengths <- unlist(lapply(1:obs.num, FUN=function(i) length(unlist(strsplit(as.list(seq.environment)[[i]],"")))))
       seq.lengths.tab <- table(seq.lengths)
       obs.have.same.length <- length(seq.lengths.tab)==1
       if(! obs.have.same.length){# DELETING THE SEQS THAT ARE DIFFERENT FROM THE MAJORITY LENGTH
          mode.1.length <- as.numeric(names(seq.lengths.tab))[which(seq.lengths.tab == max(seq.lengths.tab))[1]]
          seqs.2b.rm <- names(unlist(as.list(seq.environment)))[which(seq.lengths != mode.1.length)]
          rm(list = seqs.2b.rm ,envir = seq.environment)
          warning(paste("Observed sequences had different lengths, all but sequences of lengths", mode.1.length, "were deleted"))
       }
       seq.char <- unlist(as.list(seq.environment))
       seq.names <- names(seq.char)
       all.seq <- unlist(strsplit(paste(seq.char,sep="",collapse=""),split=""))
  
       nrow <- length(seq.char)
       ncol <- length(all.seq)/nrow
       seq.mtx <- aperm(array(all.seq,dim =c(ncol,nrow)), perm = c(2,1))
       rownames(seq.mtx) <- seq.names
       return(seq.mtx)
    }

library("seqinr")

#' Internal function getSeqMtx
#' 
#' Reads an alignment produced by \code{seqinr}'s \code{read.alignment} function and produces a \eqn{seqs x seqlen}
#' matrix of characters.  
#' 
#' The length of the first sequence in the alignment is taken to be the "true" length. If any sequence in the alignment
#' has a length that differs from the "true" length, by default, \code{getSeqMtx} throws an error.  If \code{warn.and.recover}
#' is TRUE, the offending sequence(s) will be removed from the output matrix, and a warning will be generated.  In 
#' either case, the warning or error message will contain the names of the offending sequences.
#' 
#' @param fname Name of aligned fasta file
#' @param prefix path to prefix to file name.  Defaults to <siever-package-location>/extdata
#' @param warn.and.recover Boolean.  False by default.
#' @return Matrix of characters.  Rows = number of sequences of the "correct" length.  Columns = alignment positions. Row labels = sequence names.
#' @seealso \code{\link{convert.env.to.matrix}}
getSeqMtx <- 
   function(
      fname,
      prefix=paste(system.file(package="siever"),"/extdata/",sep=""),
      warn.and.recover = FALSE
   ) {
         aln <- read.alignment(file=paste(prefix,fname,sep=""),format="fasta",forceToLower=FALSE)
         # Just the sequence part
         alnSeqs <- strsplit(aln$seq,"")
         ### TODO: Die by default if it's not an aligned fasta file.
         # Vector of indicators telling which sequence lengths 
         #   are equal to the length of first sequence
         alnGoodSeqs <- unlist(lapply(alnSeqs,function(l)length(l)==length(alnSeqs[[1]])))
         # Are some of the sequences a different length than the first sequence?
         if(!all(alnGoodSeqs)) {
            baddies <- paste(paste(aln$nam[!alnGoodSeqs],collapse=", "),"have incorrect sequence lengths")
            if(warn.and.recover) {
               warning(baddies)
            } else { 
               stop(baddies)
            }
         }      
         # Just the "good" sequencs 
         alnSeqs <- alnSeqs[alnGoodSeqs]
         # Matrix of "good" sequences; original sequences names become row labels
         return(
            matrix(unlist(alnSeqs),
                   byrow=TRUE,
                   nrow=length(alnSeqs),
                   ncol=length(alnSeqs[[1]]),
                   dimnames=list(aln$nam[alnGoodSeqs])
            )
         )
   }

#' Internal function filterAcceptable
#' 
#' Reads a matrix of characters.  Converts any character \emph{not} in \code{acceptable.chars} into \code{replace.with}.
#' 
#' @param mtx A matrix or array of single characters.
#' @param acceptable.chars Vector of characters. Defaults to legal IUPAC single-letter codes for amino acids + "-".
#' @param replace.with Single character.  Any unacceptable character is replaced with this.  Default is NA
#' 
#' @return Matrix of characters with replacements
#' @seealso \code{\link{convert.env.to.matrix}}
filterAcceptable <-
   function(
      mtx,
      acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ),
      replace.with = NA
   ) {
      unacceptable.pattern <- paste("[^",paste(acceptable.chars,collapse=""),"]",sep="")
      return(gsub(unacceptable.pattern,replace.with,mtx))
   } 

#' Internal function
#' This funciton update the site.sets.list based on the sites that will be kept according to include.site
#' @param include.site <to be described>
#' @param site.sets.list <to be described>

update.site.sets.list <- function(include.site, site.sets.list){
  
  if(!is.null(site.sets.list)){
    renumbered.sites <- cumsum(include.site) * include.site
    tmp.site.sets.list <- lapply(site.sets.list, FUN= function(.set){
      tmp.set <- renumbered.sites[.set]
      return(tmp.set[-which(tmp.set == 0)] )
    })
    site.sets.list <- tmp.site.sets.list
    rm(tmp.site.sets.list)
  }
  return(site.sets.list)
}

#' Internal function
#' This function given a vector of stats and the permuted null matrix (rows = permutations, columns = sites), it computes a vector of per site p-values
#' @param stats is the test stat
#' @param perm.stats, is the permuted null

get.perm.p.value <- function(stats, perm.stats){
  stats.comp.mtx <- matrix(stats,nrow=nrow(perm.stats),ncol = length(stats),byrow=T)
  
  non.na.count <- apply(!is.na(perm.stats),MARGIN=2,FUN=sum)
  tail.L.count <- apply(stats.comp.mtx >= perm.stats,MARGIN=2,FUN=function(.site) sum(.site,na.rm=T))
  tail.R.count <- apply(stats.comp.mtx <= perm.stats,MARGIN=2,FUN=function(.site) sum(.site,na.rm=T))
  tail.p <- 2 * apply(rbind(tail.L.count,tail.R.count),MARGIN=2,FUN=min)/non.na.count
  tail.p[is.na(stats)] <- NA
  p.value <- apply(rbind(tail.p,1),MARGIN=2,FUN=function(.site) min(.site,na.rm=T))
  return(p.value)
}

#' Internal function
#' This function computes per site p-values vector, when permuted null is not present it makes distributional assumption (e.g. t-distribution)
#' @param stats is the test stat
#' @param perm.stats, is the permuted null
#' @param type is the statistics type
get.p.value <- function(stats, perm.stats = NULL, type = c("tStat","mimic.smmb","KL","SKL")){
  
  if(type == "tStat" & is.null(perm.stats)){
      p.value <- 2*pt(-abs(stats$test.stats), df = stats$df)
      return(p.value)
  } 
  p.value <- get.perm.p.value(stats=stats$test.stats,perm.stats=perm.stats)
  return(p.value)
  
}

#' Internal function
#' this function checks if a number is a natural number
#' @param num numeric value
is.natural <- function(num){
  if(!is.numeric(num)) stop("only numeric values")
  return(1<=num & num - round(num)==0)
}