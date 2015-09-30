makeSieverTest <- function ( statFunction = SieverStatFunction.T, needsReferenceSeq = TRUE, substitutionMatrix = "default", oneSeqPerSubject = TRUE ),

SieverParameters.TESTS <- list(
    GWJ = makeSieverTest( statFunction = SieverStatFunction.T, needsReferenceSeq = TRUE, substitutionMatrix = "default", sequenceSetAggregatorFunction = NA ),
    EGWJ = makeSieverTest( statFunction = SieverStatFunction.T, needsReferenceSeq = TRUE, substitutionMatrix = "default", sequenceSetAggregatorFunction = mean ),
    KL = makeSieverTest( statFunction = SieverStatFunction.KL, needsReferenceSeq = FALSE, substitutionMatrix = NA, sequenceSetAggregatorFunction = NA ),
    SKL = makeSieverTest( statFunction = SieverStatFunction.SKL, needsReferenceSeq = FALSE, substitutionMatrix = NA, sequenceSetAggregatorFunction = NA )
    SMMB = makeSieverTest( statFunction = SieverStatFunction.SMMB, needsReferenceSeq = TRUE, substitutionMatrix = "Hamming", sequenceSetAggregatorFunction = sum )
);
SieverParameters.TEST.TYPES <- names( SieverParameters.TESTS );
### Idea: need a type-specific validateMe function that is called (only) right before tests actually run.  for now called at creation time, which is ok.

#' SieverParameters.AA.ACCEPTABLE.CHARS definition
SieverParameters.AA.ACCEPTABLE.CHARS <- c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" );

#' SieverParameters CLASS CONSTRUCTOR
#'
#' @param weightMtx A 2-dimensional square "substitution" matrix with row and column names given by single characters found in the sequence files: AAs or nucleotides.  For reference-sequence-based methods, the reference AA at each site determines the row of the matrix from which weights are determined (for each breakthrough sequence, its AA determines the column of weightMtx to use).  Leave as NULL to use the identity matrix labeled with acceptableChars, which results in Hamming distance comparisons.  Note that gap characters should be explicitly included in the matrix (with label "-"), and that any observed character not found in the matrix labels will be treated as missing data.  See also acceptableChars, which is used only if weightMtx is NULL.
#' @param acceptableChars an array of characters used when weightMtx is given as NULL to determine the length and label names of the identity weight matrix; defaults to the 20 non-ambiguous amino acid IUPAC symbols and the gap char ("-") (21 chars total).  Note that the gap character(s) should be explicitly included in this array (with label "-"), and that any observed character not found in the matrix labels will be treated as missing data.
#' @param testType One of the character strings in SieverParameters.TEST.TYPES; by default the GWJ nonparametric T test method is used, which is the statistic Z_2^A(i) given in equation 4 of Gilbert, Wu, Jobes (2008)).
#' @param siteSets A list of arrays containing character site names or numerical site indices (1-indexed) for grouping sites together (see withinSiteSetWeightAccumulationFn).  If left at the default value of NULL, each site will be its own set of cardinality one.
#' @param obs.vSeq.sets.list <to be described>
#' @param obs.pSeq.sets.list <to be described>
#' @param withinSiteSetWeightAccumulationFn <to be described>
#' @param within.obs.acc.fn <to be described>
#'
#' @export
SieverParameters <- function (
    weightMtx = NULL,
    testType = SieverParameters.TEST.TYPES,
    siteSets = NULL,
    vSeqSets = NULL,
    pSeqSets = NULL,
    acceptableChars = SieverParameters.AA.ACCEPTABLE.CHARS,
    within.site.acc.fn = sum,
    within.Seqset.acc.fn = mean,
    between.site.acc.fns = list()
) {
    
    testType <-
        tryCatch( match.arg( testType, SieverParameters.TEST.TYPES, several.ok = FALSE ), error = function( e ) {
      message( "testType needs to be in (", paste( SieverParameters.TEST.TYPES, collapse = ", " ), ")" )
      return(NA)
    },finally={})
    if(is.na(testType)){
      return()
    }
    mimic.smmb <- testType == "mimic.smmb"
    KL <- testType == "KL"
    SKL <- testType == "SKL"
    if(KL || SKL){
      siteSets <- NULL
      obs.vSeq.sets.list <- NULL
      obs.pSeq.sets.list <- NULL
      message("NOTE: KL and SKL option currently does not support accumulation within sites or observation sets")
    }

    # In case the user specified a
    if( is.null( siteSets ) ) {
        site.sets.have.cardinality.one <- TRUE;
        ## TODO: HERE replace siteSets with one-to-one map?
    } else if( !is.null( siteSets ) ) {
        site.sets.have.cardinality.one <-
            all( unlist( lapply( siteSets, FUN=function( set ) { length( set ) == 1 } ) ) );
    }
    
   if( is.null( acceptableChars ) ) {
     acceptableChars <- SieverParameters.AA.ACCEPTABLE.CHARS;
    }
    if( !is.null( weightMtx ) ) {
      acceptableChars <- setdiff( rownames( weightMtx ), 'X' )
      message("NOTE: When weight matrix is provided acceptable characters are overridden based on row names of the matrix\n") 
    }
    
    if( mimic.smmb ) {
      #within.obs.acc.fn <- function(nums) return(mean(nums, na.rm=T))
      within.obs.acc.fn <- sum
      message("NOTE: When mimic option is chosen function for accumulation-across-observations-in-a-group is always overridden to sum\n") 
    }
    
    opts <-
      list(
        weightMtx = weightMtx,
        testType = testType,
        siteSets = siteSets,
        site.sets.have.cardinality.one = site.sets.have.cardinality.one,
        obs.vSeq.sets.list = obs.vSeq.sets.list,
        obs.pSeq.sets.list = obs.pSeq.sets.list,
        acceptableChars = acceptableChars,
        withinSiteSetWeightAccumulationFn = withinSiteSetWeightAccumulationFn,
        within.obs.acc.fn = within.obs.acc.fn,
        call = match.call()
      )
    
    class(opts) <- "SieverParameters"
    return(opts)
}  

#### TODO: UPDATE

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
