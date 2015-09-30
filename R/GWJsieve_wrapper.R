##' THIS IS A WRAPPER FOR \code{\link{GWJsieve}} TO PROVIDE THE SAME INTERFACE AS (the obsolescent) \code{tStat}.
#'
#'  NOTE: NOT ALL OPTIONS ARE IMPLEMENTED AT THIS POINT
#'  Many of the parameters in this shim are shipped to \code{\link{GWJsieve}} via the \code{\link{weightopt}} 
#'  object in the current sieveR package. 
#' @param vaccine.seq.envir Environment containing aligned protein sequences for vaccine recipients
#' @param placebo.seq.envir Environment containing aligned protein sequences for placebo recipients
#' @param insert.seq.envir Environment containing viral insert protein sequences
#' @param weight.matrix Matrix containing amino-acid substitution scores
#' @param site.sets.list List of protein sites to analyze
#' @param use.f.test Boolean
#' @param return.t.test.result Boolean
#' @param weights.across.sites.in.a.set.init <to be described>
#' @param weights.across.sites.in.a.set.accumulation.fn <to be described>
#' @param vaccine.sequence.sets.list <to be described>
#' @param placebo.sequence.sets.list <to be described>
#' @param mimic.smmb Boolean
#' @param weights.across.sequences.in.a.set.accumulation.fn <to be described>
#' @param acceptable.chars Alphabet for representing amino acids and gaps
#' @param instead.return.weights Boolean
#' @keywords sieve-analysis T-statistic GWJ
#' @export
#'
GWJsieve.wrapper <-
  function(
    vaccine.seq.envir,
    placebo.seq.envir,
    insert.seq.envir,
    weight.matrix = NULL,
    site.sets.list = NULL,
    use.f.test = FALSE,
    return.t.test.result = FALSE,
    weights.across.sites.in.a.set.init = 0,
    weights.across.sites.in.a.set.accumulation.fn = sum,
    vaccine.sequence.sets.list = NULL,
    placebo.sequence.sets.list = NULL,
    mimic.smmb = FALSE,
    weights.across.sequences.in.a.set.accumulation.fn =
      if( mimic.smmb ) { sum } else { mean },
    acceptable.chars =
      if( is.null( weight.matrix ) ) {
        c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" )
      } else {
        setdiff( rownames( weight.matrix ), 'X' )
      },
    instead.return.weights = FALSE ) {
    
    insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]]
    vaccine.seq.mtx <- convert.env.to.matrix(seq.environment=vaccine.seq.envir)
    placebo.seq.mtx <- convert.env.to.matrix(seq.environment=placebo.seq.envir)
    stat.type <- "tStat"
    if(mimic.smmb) stat.type <- "mimic.smmb"
    
    wOpt <-
      weightopt(
        weight.mtx = weight.matrix,
        stat.type = stat.type,
        insert.char.vector,
        site.sets.list= site.sets.list,
        obs.vSeq.sets.list = vaccine.sequence.sets.list,
        obs.pSeq.sets.list = placebo.sequence.sets.list,
        acceptable.chars = acceptable.chars,
        within.sites.acc.fn = weights.across.sites.in.a.set.accumulation.fn,
        within.obs.acc.fn = weights.across.sequences.in.a.set.accumulation.fn
      )
    
    data.mtx <- rbind(vaccine.seq.mtx, placebo.seq.mtx)
    v.indices <- 1:nrow(vaccine.seq.mtx)
    out <- GWJsieve(data.mtx, v.indices, wOpt = wOpt,sw.test=T)
    return(out)
  }