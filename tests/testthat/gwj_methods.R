#################################################################################
### Description: Performs testing procedures to identify type A (t test) and type B (Kullback-Leibler test)
###              signature sites based on Gilbert, Wu, and Jobes (2008)
### Author: Michal Juraska (using a few functions from Youyi Fong and Paul Edlefsen)
### Date:   Aug 11, 2011
### Examples:
### KULLBACK-LEIBLER TEST
### 1) Compute permutation-based unadjusted per position p-values using the Kullback-Leibler test:
### 1a) All sites:
### klTestPermPval("vaccine_maj_gag_00.fasta", "placebo_maj_gag_00.fasta", num.perm=100)
###
### 1b) Sites selected by the TRUE/FALSE vector 'include.site':
### site <- logical(542); site[501:542] <- TRUE
### klTestPermPval("vaccine_maj_gag_00.fasta", "placebo_maj_gag_00.fasta", include.site=site, num.perm=100)
###
### 2) In addition to 1a), determine whether we reject the null hypothesis based on fine-tuning the significance threshold as described in Section 3.3 of GWJ (2008):
### klTestPermPval("vaccine_maj_gag_00.fasta", "placebo_maj_gag_00.fasta", num.perm=100, determine.signif=TRUE, num.perm.signif=100, num.FP.rejections=20)
###
### T TEST
### 3) Compute permutation-based unadjusted per position p-values using the t test:
### 3a) All sites:
### tTestPermPval("vaccine_maj_gag_00.fasta", "placebo_maj_gag_00.fasta", "insert_gag_00.fasta", "PAMamong.txt", num.perm=100)
###
### 3b) Sites selected by the TRUE/FALSE vector 'include.site':
### site <- logical(542); site[501:542] <- TRUE
### tTestPermPval("vaccine_maj_gag_00.fasta", "placebo_maj_gag_00.fasta", "insert_gag_00.fasta", "PAMamong.txt", include.site=site, num.perm=100)
###
##################################################################################

source( "sequence.util.R" )

######################### Insert-dependent ###############
### 'tStat' returns values of the t test statistic per aa position (or, if site.sets.list is non-NULL, one stat per site-set).
# If weight.matrix is NULL, weights will be 0 for insert-match and 1 for insert-mismatch.
# If site.sets.list is non-NULL, it should be a list of vectors of site numbers, with numbers referring to positions in the given seq.envirs. -- Note that these numbers must be in the range 1..nchar( get( ls( env=insert.seq.envir )[1], env=insert.seq.envir ) ), so if you've done any screening via eg screenOutSites(..), you'll need to make sure to update the values.  See screenOutSites(..) and updateSiteSetsList(..).
# Note that by default the weights will be summed across the sites in a set, but you can change this behavior by changing the weights.across.sites.in.a.set.accumulation.fn -- but note that you should also change weights.across.sites.in.a.set.init from 0!  This is the initial weight -- so for instance if you're using weights.across.sites.in.a.set.accumulation.fn = prod, you should set weights.across.sites.in.a.set.init to 1; if you're using weights.across.sites.in.a.set.accumulation.fn = min, you should set weights.across.sites.in.a.set.init to Inf; if you're using max, set weights.across.sites.in.a.set.init to -Inf.
#### ERE WE ARE.  Working on implementing the use.f.test option to mimic the fortran code implementing the GWJ "Euclidean" stat as seen in GWJ 2006 (and GWJ 2008, but there it's missing the weight matrix).  The use.f.test option, if TRUE, should mimic the option in that fortran code with "qq" = 2 (which happens when the Blosum weight matrix does NOT have a 0 in the 0,0 position.  It's the Z_E stat in GWJ 2006 (fmla 1), except (currently) without the lambda1 value (or equivalently with lambda1=0).  It's the Z_{E2}^B stat in GWJ 2008 except that has no weight matrix.
#### TODO: DEBUG/TEST the use.f.test option.
#### NOTE: With the use.f.test value at its default (FALSE), this implements the t-test version of the stat, which is GWJ 2008 Z_2^A in equation 4.
#### TODO: NOTE: FOR NOW, the f test option can't be used in conjunction with mimic.smmb, return.t.test.result or with site sets or sequence sets (of cardinality greater than 1) !!!!!
tStat <- function( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, weight.matrix = NULL, site.sets.list = NULL, use.f.test = FALSE, return.t.test.result = FALSE, weights.across.sites.in.a.set.init = 0, weights.across.sites.in.a.set.accumulation.fn = sum, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, mimic.smmb = FALSE, weights.across.sequences.in.a.set.accumulation.fn = if( mimic.smmb ) { sum } else { mean }, acceptable.chars = if( is.null( weight.matrix ) ) { c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" )} else{ setdiff( rownames( weight.matrix ), 'X' ) }, instead.return.weights = FALSE )
{
  insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]];

  site.sets.have.cardinality.one <- NA;
  if( is.null( site.sets.list ) ) {
    site.sets.list = lapply( 1:length( insert.char.vector ), function( .site ) { .site } );
    site.sets.have.cardinality.one <- TRUE;
  }
  num.site.sets <- length( site.sets.list );

  vaccine.sequence.sets.have.cardinality.one <- NA;
  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.seq.envir ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.seq.envir );
    vaccine.sequence.sets.have.cardinality.one <- TRUE;
  }
  placebo.sequence.sets.have.cardinality.one <- NA;
  if( is.null( placebo.sequence.sets.list ) ) {
    stopifnot( !is.null( placebo.seq.envir ) );
    placebo.sequence.sets.list = lapply( ls( env=placebo.seq.envir ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.seq.envir );
    placebo.sequence.sets.have.cardinality.one <- TRUE;
  }

  ## Check that the sites and sequence sets all have cardinality 1 for use with the use.f.test option.  Also check that return.t.test.result and mimic.smmb are FALSE.
  if( use.f.test ) {
    stopifnot( !return.t.test.result );
    stopifnot( !mimic.smmb );
    # CHECK that the sequence sets and site sets all have cardinality 1.  Otherwise die and complain bitterly.
    if( is.na( site.sets.have.cardinality.one ) ) {
      site.sets.have.cardinality.one <- all( sapply( site.sets.list, function( .element ){ length( .element ) == 1 } )  );
      stopifnot( site.sets.have.cardinality.one );
    }
    if( is.na( vaccine.sequence.sets.have.cardinality.one ) ) {
      vaccine.sequence.sets.have.cardinality.one <- all( sapply( vaccine.sequence.sets.list, function( .element ){ length( .element ) == 1 } ) );
      stopifnot( vaccine.sequence.sets.have.cardinality.one );
    }
    if( is.na( placebo.sequence.sets.have.cardinality.one ) ) {
      placebo.sequence.sets.have.cardinality.one <- all( sapply( placebo.sequence.sets.list, function( .element ){ length( .element ) == 1 } ) );
      stopifnot( placebo.sequence.sets.have.cardinality.one );
    }
  } # End if use.f.test
  
  .vaccine.weights.by.pos <- list();
  .placebo.weights.by.pos <- list();
  if( mimic.smmb ) {
    .vaccine.observed.seqs.count.by.pos <- list();
    .placebo.observed.seqs.count.by.pos <- list();
  }
  # Start the clock!
  #ptm <- proc.time()
  for( .pos.i in unique( unlist( site.sets.list ) ) ) {
    
    if( is.na( .pos.i ) ) {
      next;
    }
    insert.char.pos <- insert.char.vector[ .pos.i ];
    
    vaccine.chars.pos <- aaPerPosition( vaccine.seq.envir, .pos.i );
    vaccine.chars.pos[ !( vaccine.chars.pos %in% acceptable.chars ) ] <- NA;
    .seq.weights <- rep( NA, length( vaccine.chars.pos ) );
    names( .seq.weights ) <- names( vaccine.chars.pos );
    for( .seq.k in 1:length( .seq.weights ) ) {
      if( is.na( vaccine.chars.pos[ .seq.k ] ) ) {
        .weight <- NA;
      } else if( is.null( weight.matrix ) ) {
        # Use 0/1 indicators of insert mismatch (1 iff mismatch)
        .weight <- vaccine.chars.pos[ .seq.k ] != insert.char.pos;
      } else {
        .weight <- weight.matrix[ insert.char.pos, vaccine.chars.pos[ .seq.k ] ];
      }
      .seq.weights[ .seq.k ] <- .weight;
    } # End foreach vaccine seq.k
    .vaccine.weights.by.pos[[ .pos.i ]] <- sapply( vaccine.sequence.sets.list, function( .sequences.in.set ) { weights.across.sequences.in.a.set.accumulation.fn( .seq.weights[ .sequences.in.set ], na.rm = T ) } );
    if( mimic.smmb ) {
      .vaccine.observed.seqs.count.by.pos[[ .pos.i ]] <- sapply( vaccine.sequence.sets.list, function( .sequences.in.set ) { sum( !is.na( .seq.weights[ .sequences.in.set ] ) ) } );
    }

    if( !is.null( placebo.seq.envir ) ) {
      placebo.chars.pos <- aaPerPosition( placebo.seq.envir, .pos.i );
      placebo.chars.pos[ !( placebo.chars.pos %in% acceptable.chars ) ] <- NA;
      .seq.weights <- rep( NA, length( placebo.chars.pos ) );
      names( .seq.weights ) <- names( placebo.chars.pos );
      for( .seq.k in 1:length( .seq.weights ) ) {
        if( is.na( placebo.chars.pos[ .seq.k ] ) ) {
          .weight <- NA;
        } else if( is.null( weight.matrix ) ) {
          # Use 0/1 indicators of insert mismatch (1 iff mismatch)
          .weight <- as.numeric( placebo.chars.pos[ .seq.k ] != insert.char.pos );
        } else {
          .weight <- weight.matrix[ insert.char.pos, placebo.chars.pos[ .seq.k ] ];
        }
        .seq.weights[ .seq.k ] <- .weight;
      } # End foreach placebo seq.k
    } # End if !is.null( placebo.seq.envir )
    .placebo.weights.by.pos[[ .pos.i ]] <- sapply( placebo.sequence.sets.list, function( .sequences.in.set ) { weights.across.sequences.in.a.set.accumulation.fn( .seq.weights[ .sequences.in.set ], na.rm = T ) } );
    if( mimic.smmb ) {
      .placebo.observed.seqs.count.by.pos[[ .pos.i ]] <- sapply( placebo.sequence.sets.list, function( .sequences.in.set ) { sum( !is.na( .seq.weights[ .sequences.in.set ] ) ) } );
    }
  } # End foreach .pos.i, calc weights..
  # Stop the clock
  #print(proc.time() - ptm)
  
  test.statistics <- rep( NA, num.site.sets );
  if( instead.return.weights ) {
    test.weights <- list();
  }
  
  for( site.set.j in 1:num.site.sets ) {
    if( all( is.na( site.sets.list[[ site.set.j ]] ) ) ) {
      # No sites in this set!
      next;
    }
    vaccine.weights <- rep( weights.across.sites.in.a.set.init, length( vaccine.sequence.sets.list ) );
    placebo.weights <- rep( weights.across.sites.in.a.set.init, length( placebo.sequence.sets.list ) );
    num.vaccine.weights <- sum( !is.na( vaccine.weights ) );
    num.placebo.weights <- sum( !is.na( placebo.weights ) );
    num.weights <- num.vaccine.weights + num.placebo.weights;
    for( .pos.i in site.sets.list[[ site.set.j ]] ) {
      if( is.na( .pos.i ) ) {
        # Skip NAs.
        next;
      }
      vaccine.weights <- apply( cbind( vaccine.weights, .vaccine.weights.by.pos[[ .pos.i ]] ), 1, weights.across.sites.in.a.set.accumulation.fn );
      placebo.weights <- apply( cbind( placebo.weights, .placebo.weights.by.pos[[ .pos.i ]] ), 1, weights.across.sites.in.a.set.accumulation.fn );
    } # End foreach pos .pos.i in site.set.j, update weights.
    if( instead.return.weights ) {
      test.weights[[ site.set.j ]] <- list( vaccine.weights = vaccine.weights, placebo.weights = placebo.weights );
    }
    
    # Special cases: if there's no variation in either group, the test stat should revert to the difference between the group constant-values (not Studentized means) - but if there's no variation _anywhere_ (all values equal a constant) then return NA to indicate that there's no point in doing permutations or bootstraps.
    ## TODO: REMOVE. DEBUGGING.
    if( is.na( all( any( all( is.na( vaccine.weights ) ), all( vaccine.weights == vaccine.weights[ 1 ] ), na.rm = T ), any( all( is.na( placebo.weights ) ), all( placebo.weights == placebo.weights[ 1 ] ), na.rm = T ) ) ) ) {
      print( "NA: all( any( all( is.na( vaccine.weights ) ), all( vaccine.weights == vaccine.weights[ 1 ] ), na.rm = T ), any( all( is.na( placebo.weights ) ), all( placebo.weights == placebo.weights[ 1 ] ), na.rm = T ) )" );
      print( all( is.na( vaccine.weights ) ) );
      print( all( vaccine.weights == vaccine.weights[ 1 ] ) );
      print( all( is.na( placebo.weights ) ) );
      print( all( placebo.weights == placebo.weights[ 1 ] ) );
      stop();
    }
    if( all( any( all( is.na( vaccine.weights ) ), all( vaccine.weights == vaccine.weights[ 1 ] ), na.rm = T ), any( all( is.na( placebo.weights ) ), all( placebo.weights == placebo.weights[ 1 ] ), na.rm = T ) ) ) {
      if( ( is.na( placebo.weights[ 1 ] ) || is.na( vaccine.weights[ 1 ] ) ) ||
          ( placebo.weights[ 1 ] == vaccine.weights[ 1 ] ) ) {
        test.statistics[ site.set.j ] <- NA;
      } else {
        test.statistics[ site.set.j ] <-
          mean( vaccine.weights, na.rm = T ) - mean( placebo.weights, na.rm = T );
      }
    } else if( use.f.test ) {
      ## TODO: Allow use of (asymptotic, sometimes) F distribution?  For now we force using permutation.
      test.statistics[ site.set.j ] <- 0; # We'll be summing, so we start at 0
      .num.nonzero.components <- 0;
      vaccine.chars.pos <- aaPerPosition( vaccine.seq.envir, site.set.j );
      placebo.chars.pos <- aaPerPosition( placebo.seq.envir, site.set.j );
      insert.char.pos <- insert.char.vector[ site.set.j ]; 
      for( .category.a in 1:length( acceptable.chars ) ) {
        ### TODO: Implement use of lambda1 (see GWJ 2006 paper, fmla 1)?  Note that in the genomescan.DataAnalysis.f fortran code implementing this for the Step analysis, the lambda1 value was not used (or it was 0, depending on how you think of it).  This seems to have been decided based on the simulation results from genomescan.simulate.11.05 (OR genomescansimulate.8.2006.f?).
        .vaccine.vector.of.indicators.of.having.category.a <- ( vaccine.chars.pos == acceptable.chars[ .category.a ] );
        .placebo.vector.of.indicators.of.having.category.a <- ( placebo.chars.pos == acceptable.chars[ .category.a ] );
        .vaccine.prob.of.category.a <- mean( .vaccine.vector.of.indicators.of.having.category.a, na.rm=T );
        .placebo.prob.of.category.a <- mean( .placebo.vector.of.indicators.of.having.category.a, na.rm=T );
        .v <- weight.matrix[ insert.char.pos, acceptable.chars[ .category.a ] ] ** 2;
        .v <- .v * (num.vaccine.weights - 1)* .vaccine.prob.of.category.a *(1 - .vaccine.prob.of.category.a)/(num.weights - 2) + (num.placebo.weights - 1)* .placebo.prob.of.category.a *(1 - .placebo.prob.of.category.a)/(num.weights - 2);
        if( .v > 0 ) {
          test.statistics[ site.set.j ] <- test.statistics[ site.set.j ] + ( ( mean( vaccine.weights * .vaccine.vector.of.indicators.of.having.category.a ) - mean( placebo.weights * .placebo.vector.of.indicators.of.having.category.a ) ) ** 2 ) / .v;
          .num.nonzero.components <- .num.nonzero.components + 1;
        }
      } # End foreach .category.a ..
      if( .num.nonzero.components != 0 ) {
        # Gotta get the constant (called C_E in the 2006 version of GWJ) in there..
        .constant <- ( ( num.weights - 2 - .num.nonzero.components + 1 ) / ( .num.nonzero.components * ( num.weights - 2 ) ) );
        .constant <- .constant * ( ( num.vaccine.weights - 1 ) * ( num.placebo.weights - 1 ) / ( num.weights - 2 ) );
        test.statistics[ site.set.j ] <- test.statistics[ site.set.j ] * .constant;
      } else {
        stopifnot( test.statistics[ site.set.j ] == 0 );
      }
    } else { # if use.f.test .. else ..
      if( return.t.test.result ) {
        t.test.result <- t.test( vaccine.weights, placebo.weights );
        test.statistics[ site.set.j ] <- list( t.test.result = t.test.result );
      } else {
        if( mimic.smmb ) {
          vaccine.observed.seqs.count <- sum( unlist( .vaccine.observed.seqs.count.by.pos[ site.sets.list[[ site.set.j ]] ] ) );
          if( vaccine.observed.seqs.count > 0 ) {
            phat.v <- ( sum( vaccine.weights, na.rm=T ) / vaccine.observed.seqs.count );
          } else {
            phat.v <- 0;
          }
          placebo.observed.seqs.count <- sum( unlist( .placebo.observed.seqs.count.by.pos[ site.sets.list[[ site.set.j ]] ] ) );
          if( placebo.observed.seqs.count > 0 ) {
            phat.p <- ( sum( placebo.weights, na.rm=T ) / placebo.observed.seqs.count );
          } else {
            phat.p <- 0;
          }
          test.statistics[ site.set.j ] <- ( phat.v - phat.p );
        } else {
          se <- sqrt( ((num.vaccine.weights - 1)*var(vaccine.weights, na.rm=T)/(num.weights - 2) + (num.placebo.weights - 1)*var(placebo.weights, na.rm=T)/(num.weights - 2)) * (1/num.vaccine.weights + 1/num.placebo.weights) );
          test.statistics[ site.set.j ] <- ( ( mean( vaccine.weights, na.rm=T ) - mean( placebo.weights, na.rm=T ) ) / se );
        }
      }
    } # End if use.f.test .. else ..
  } # End foreach site.set.j

  if( instead.return.weights ) {
    return( test.weights );
  }
  return( test.statistics );
} ## tStat (..)

### 'tPermTestStatistics' returns a matrix of t test statistics per aa position from 'num.perm' permuted data sets
# If weight.matrix is NULL, weights will be 0 for insert-match and 1 for insert-mismatch.
tPermTestStatistics <- function( vaccine.seq.env, placebo.seq.env, insert.seq.env, weight.matrix = NULL, site.sets.list = NULL, use.f.test = FALSE, weights.across.sites.in.a.set.init = 0, weights.across.sites.in.a.set.accumulation.fn = sum, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, mimic.smmb = FALSE, weights.across.sequences.in.a.set.accumulation.fn = if( mimic.smmb ) { sum } else { mean }, num.perm = 100, use.bootstrap = ifelse( mimic.smmb, TRUE, FALSE ), be.verbose = FALSE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{
  if( is.null( site.sets.list ) ) {
    insert.char.vector <- strsplit(as.list(insert.seq.env)[[1]],"")[[1]];
    site.sets.list = lapply( 1:length( insert.char.vector ), function( .site ) { .site } );
  }
  num.site.sets <- length( site.sets.list );

  combined.seq.env <- new.env( hash = TRUE );
  lapply( ls( env=vaccine.seq.env ), function( .vaccine.sequence.name ) { assign( .vaccine.sequence.name, get( .vaccine.sequence.name, env=vaccine.seq.env ), env=combined.seq.env ) } );
  lapply( ls( env=placebo.seq.env ), function( .placebo.sequence.name ) { assign( .placebo.sequence.name, get( .placebo.sequence.name, env=placebo.seq.env ), env=combined.seq.env ) } );
  
  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.seq.env ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.seq.env );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    placebo.sequence.sets.list = lapply( ls( env=placebo.seq.env ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.seq.env );
  }

  .renamed.vaccine.sequence.sets.list <- vaccine.sequence.sets.list;
  names( .renamed.vaccine.sequence.sets.list ) <- paste( "vaccine:", names( vaccine.sequence.sets.list ), sep = "" );
  .renamed.placebo.sequence.sets.list <- placebo.sequence.sets.list;
  names( .renamed.placebo.sequence.sets.list ) <- paste( "placebo:", names( placebo.sequence.sets.list ), sep = "" );
  
  all.sequence.sets.list <- c( .renamed.vaccine.sequence.sets.list, .renamed.placebo.sequence.sets.list );
  
  vaccine.num.sequence.sets <- length( vaccine.sequence.sets.list );
  placebo.num.sequence.sets <- length( placebo.sequence.sets.list );
  num.sequence.sets <- vaccine.num.sequence.sets + placebo.num.sequence.sets;
  stopifnot( num.sequence.sets == length( all.sequence.sets.list ) );

  test.statistics.perm <- matrix( 0, nrow=num.site.sets, ncol=num.perm );
  for( perm.i in 1:num.perm ) {
    if( be.verbose && ( ( num.perm <= 100 ) || ( ( perm.i %% floor( num.perm / 100 ) ) == 0 ) ) ) {
      cat( paste( "Calling tStat(..) for", ifelse( use.bootstrap, "bootstrap", "permutation" ), perm.i, "of", num.perm, "..\n" ) );
    }
    .permutation <- sample( 1:num.sequence.sets, replace = use.bootstrap );

    .perm.vaccine.sequence.sets.list <- all.sequence.sets.list[ .permutation[ 1:vaccine.num.sequence.sets ] ];
    .perm.placebo.sequence.sets.list <- all.sequence.sets.list[ .permutation[ ( vaccine.num.sequence.sets + 1 ):num.sequence.sets ] ];
    test.statistics.perm[ , perm.i ] <- tStat( combined.seq.env, NULL, insert.seq.env, weight.matrix=weight.matrix, site.sets.list=site.sets.list, use.f.test=use.f.test, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, vaccine.sequence.sets.list=.perm.vaccine.sequence.sets.list, placebo.sequence.sets.list=.perm.placebo.sequence.sets.list, mimic.smmb=mimic.smmb, weights.across.sequences.in.a.set.accumulation.fn=weights.across.sequences.in.a.set.accumulation.fn, acceptable.chars=acceptable.chars );
    # test.statistics.perm[ , perm.i ] <- tStat( combined.seq.env, NULL, insert.seq.env, weight.matrix=weight.matrix, site.sets.list=site.sets.list, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, vaccine.sequence.sets.list=.perm.vaccine.sequence.sets.list, placebo.sequence.sets.list=.perm.placebo.sequence.sets.list, mimic.smmb=mimic.smmb, weights.across.sequences.in.a.set.accumulation.fn=weights.across.sequences.in.a.set.accumulation.fn, acceptable.chars=acceptable.chars );
  } # End foreach permutation perm.i
  return( test.statistics.perm );
} # tPermTestStatistics (..)

### 'tTestPermPval' returns a list with the following components:
### 'p.values'               - numeric vector of unadjusted permutation-based t test p-values per aa position
# NOTE if include.site is "variable", then chooseSitesToScreen.minimumVariability(..) will be used (using the default values of min 3 insert-char and 3 non-insert-char seqs).
# NOTE: you can supply environments in place of filenames if you have already read in the fasta file(s).
# NOTE: If weight.matrix is NULL, weights will be 0 for insert-match and 1 for insert-mismatch.
# NOTE: by default each site is evaluated separately; change this behavior by specifying site sets (NULL means each site is in its own set). By default the weights will be summed across the sites in a set, but you can change this behavior by changing the weights.across.sites.in.a.set.accumulation.fn -- but note that you should also change weights.across.sites.in.a.set.init from 0!  This is the initial weight -- so for instance if you're using weights.across.sites.in.a.set.accumulation.fn = prod, you should set weights.across.sites.in.a.set.init to 1; if you're using weights.across.sites.in.a.set.accumulation.fn = min, you should set weights.across.sites.in.a.set.init to Inf; if you're using max, set weights.across.sites.in.a.set.init to -Inf.
# NOTE: by default each sequence contributes a weight individually; change this behavior for eg multiple sequences-per-subject by specifying sequence sets (NULL means each sequence is in its own set). Note that permutation to determine the null will be by set (ie by subject), not by individual sequence.  By default the weights will be averaged across the sequence in a set, but you can change this behavior by changing the weights.across.sequences.in.a.set.accumulation.fn, which should take a vector and should accept the "na.rm" argument (at least it shouldn't mind it; it'll always be called with na.rm = TRUE).
tTestPermPval <- function ( vaccine.fasta.file.name, placebo.fasta.file.name, insert.fasta.file.name, weight.matrix.file.name = NULL, weight.matrix.is.probability.matrix = FALSE, include.site = NULL, site.sets.list = NULL, use.f.test = FALSE, weights.across.sites.in.a.set.init = 0, weights.across.sites.in.a.set.accumulation.fn = sum, stats.across.site.sets.accumulation.fns = list(), vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, mimic.smmb = FALSE, weights.across.sequences.in.a.set.accumulation.fn = if( mimic.smmb ) { sum } else { mean }, num.perm = 1000, use.t.distribution.instead.of.permutation=( num.perm == 0 ), use.bootstrap = ifelse( mimic.smmb, TRUE, FALSE ), be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{
  if( is.null( weight.matrix.file.name ) ) {
    # If weight.matrix is NULL, weights will be 0 for insert-match and 1 for insert-mismatch.
    weight.matrix <- NULL;
    if( be.verbose ) {
      cat( "Using mismatch indicators instead of weights.", fill = T );
    }
  } else {
    if( weight.matrix.is.probability.matrix ) {
      # In probability form
      if( is.matrix( weight.matrix.file.name ) ) {
        weight.matrix <- 1 - weight.matrix.file.name;   ## convert to a dissimilarity matrix
        if( be.verbose ) {
          cat( "Using dissimilarity matrix converted from the provided probability-form substitution similarity matrix.", fill = T );
        }
      } else {
        weight.matrix <- 1 - readBetty(weight.matrix.file.name);   ## convert to a dissimilarity matrix
        if( be.verbose ) {
          cat( "Using dissimilarity matrix converted from the probability-form substitution similarity matrix in file:", weight.matrix.file.name, fill = T );
        }
      }
    } else {
      # In log-likelihood-ratio (etc) form
      if( is.matrix( weight.matrix.file.name ) ) {
        weight.matrix <- 0 - weight.matrix.file.name;    ## convert to a dissimilarity matrix
        if( be.verbose ) {
          cat( "Using dissimilarity matrix converted from the provided substitution similarity matrix.", fill = T );
        }
      } else {
        weight.matrix <- 0 - readPAM(weight.matrix.file.name);    ## convert to a dissimilarity matrix
        if( be.verbose ) {
          cat( "Using dissimilarity matrix converted from the substitution similarity matrix in file:", weight.matrix.file.name, fill = T );
        }
      }
    }
  } # End if is.null( weight.matrix.file.name ) .. else ..
  
  if( is.environment( insert.fasta.file.name ) ) {
    insert.env <- insert.fasta.file.name;
  } else {
    insert.env <- readFastaFile( insert.fasta.file.name, be.verbose=be.verbose, .acceptable.characters = acceptable.chars );
  }
  if( is.environment( vaccine.fasta.file.name ) ) {
    vaccine.env <- vaccine.fasta.file.name;
  } else {
    vaccine.env <- readFastaFile( vaccine.fasta.file.name, be.verbose=be.verbose, .acceptable.characters = acceptable.chars );
  }
  if( is.environment( placebo.fasta.file.name ) ) {
    placebo.env <- placebo.fasta.file.name;
  } else {
    placebo.env <- readFastaFile( placebo.fasta.file.name, be.verbose=be.verbose, .acceptable.characters = acceptable.chars );
  }

  ## TODO: ? .. For now, if include.site is _any_ character string, then use the minimumVariabilty rules.
  if( is.character( include.site ) ) {
    if( be.verbose ) {
      cat( "Calling chooseSitesToScreen.minimumVariability(..)", fill = T );
    }
    include.site <- chooseSitesToScreen.minimumVariability( vaccine.env, placebo.env, insert.env, vaccine.sequence.sets.list=vaccine.sequence.sets.list, placebo.sequence.sets.list=placebo.sequence.sets.list, acceptable.chars=acceptable.chars );
  }
  
  if( !is.null( include.site ) ) {
    insert.env <- screenOutSites( insert.env, include.site );
    vaccine.env <- screenOutSites( vaccine.env, include.site );
    placebo.env <- screenOutSites( placebo.env, include.site );
    if( is.null( site.sets.list ) ) {
      site.sets.list <- lapply( which( include.site ), function( .site ) { .site } );
      names( site.sets.list ) <- which( include.site );
    }
    site.sets.list <- updateSiteSetsList( include.site, site.sets.list );
  }

  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.env ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.env );
  } else {
    # I haven't yet implemented the use.t.distribution.instead.of.permutation option for the case where there's multiple seqs per subject.
    stopifnot( !use.t.distribution.instead.of.permutation );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    stopifnot( !is.null( placebo.env ) );
    placebo.sequence.sets.list = lapply( ls( env=placebo.env ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.env );
  } else {
    # I haven't yet implemented the use.t.distribution.instead.of.permutation option for the case where there's multiple seqs per subject.
    stopifnot( !use.t.distribution.instead.of.permutation );
  }
  
  ### NOTE that these are two-sided tests!
  if( use.t.distribution.instead.of.permutation ) {
    pvals <- sapply( tStat( vaccine.env, placebo.env, insert.env, weight.matrix, site.sets.list, return.t.test.result=TRUE, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, acceptable.chars=acceptable.chars ), function( .t.test.result ) { if( all( is.na( .t.test.result ) ) ) { return( NA ); } else { return( .t.test.result$p.value ); } } );
  } else {
    if( be.verbose ) {
      cat( "Calling tStat(..) for main test stats..\n" );
    }
    test.statistics <- tStat( vaccine.env, placebo.env, insert.env, weight.matrix, site.sets.list=site.sets.list, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, vaccine.sequence.sets.list=vaccine.sequence.sets.list, placebo.sequence.sets.list=placebo.sequence.sets.list, mimic.smmb=mimic.smmb, weights.across.sequences.in.a.set.accumulation.fn=weights.across.sequences.in.a.set.accumulation.fn, acceptable.chars=acceptable.chars );
    if( !is.null( stats.across.site.sets.accumulation.fns ) && ( length( stats.across.site.sets.accumulation.fns ) > 0 ) ) {
      global.stats <- lapply( stats.across.site.sets.accumulation.fns, function( .fn ) { .fn( test.statistics ) } );
    } else {
      global.stats <- NA;
    }
    if( be.verbose ) {
      cat( "Calling tStat(..) for", ifelse( use.bootstrap, "bootstraps", "permutations" ), "test stats..\n" );
    }
    test.statistics.perm <- tPermTestStatistics( vaccine.env, placebo.env, insert.env, weight.matrix, num.perm=num.perm, site.sets.list=site.sets.list, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, vaccine.sequence.sets.list=vaccine.sequence.sets.list, placebo.sequence.sets.list=placebo.sequence.sets.list, mimic.smmb=mimic.smmb, weights.across.sequences.in.a.set.accumulation.fn=weights.across.sequences.in.a.set.accumulation.fn, use.bootstrap=use.bootstrap, be.verbose=be.verbose, acceptable.chars=acceptable.chars );
    ##### NOTE: The tTest stat should be symmetric about 0 under the null but I've not verified it.  Thus we do not use this:
    #pvals <- apply( abs(test.statistics.perm) >= abs(test.statistics), 1, mean );
    ##### NOTE: instead we use this:
    calculatePValue <- function ( .pos.i, add.one.for.observed.test.statistic = FALSE ) { if( is.na( test.statistics[ .pos.i ] ) ) { return( NA ); } else { return( min( 1, 2 * ( min( sum( test.statistics[ .pos.i ] >= test.statistics.perm[ .pos.i,  ], na.rm=T ), sum( test.statistics[ .pos.i ] <= test.statistics.perm[ .pos.i,  ], na.rm=T ) ) + as.numeric( add.one.for.observed.test.statistic ) ) / (  sum( !is.na( test.statistics.perm[ .pos.i,  ] ) ) + as.numeric( add.one.for.observed.test.statistic ) ), na.rm=T ) ); } };
    pvals <- sapply( 1:nrow( test.statistics.perm ), calculatePValue );
    ## TODO: REMOVE
    if( any( is.na( pvals ) ) ) {
      # This could happen if the weights end up having no variation, which could happen even when the AAs do have variation (as usually assured by our minimum variation screen), eg when the weight.matrix contains indicators of a categorical property of the AAs.
      print( "Warning: Got NA pval for the following sites." )
      print( which( is.na( pvals ) ) );
      # na.pval.pos <- which( is.na( pvals ) )[ 1 ];
      # print( test.statistics[ na.pval.pos ] );
      # ## Recalculate it but get the weights.
      # test.weights <- tStat( vaccine.env, placebo.env, insert.env, weight.matrix, site.sets.list=site.sets.list, weights.across.sites.in.a.set.init=weights.across.sites.in.a.set.init, weights.across.sites.in.a.set.accumulation.fn=weights.across.sites.in.a.set.accumulation.fn, vaccine.sequence.sets.list=vaccine.sequence.sets.list, placebo.sequence.sets.list=placebo.sequence.sets.list, mimic.smmb=mimic.smmb, weights.across.sequences.in.a.set.accumulation.fn=weights.across.sequences.in.a.set.accumulation.fn, acceptable.chars=acceptable.chars, instead.return.weights = TRUE );
      # print( test.weights[[ na.pval.pos ]] );
      # #print( test.statistics.perm[ na.pval.pos, ] );
    }
    calculatePValue.generic <- function ( .test.statistic, .null.test.statistics, add.one.for.observed.test.statistic = FALSE ) { if( is.na( .test.statistic ) ) { NA } else { min( 1, 2 * ( min( sum( .test.statistic >= .null.test.statistics, na.rm=T ), sum( .test.statistic <= .null.test.statistics, na.rm=T ) ) + as.numeric( add.one.for.observed.test.statistic ) ) / (  sum( !is.na( .null.test.statistics ) ) + as.numeric( add.one.for.observed.test.statistic ) ), na.rm=T ) } };
    if( !is.null( stats.across.site.sets.accumulation.fns ) && ( length( stats.across.site.sets.accumulation.fns ) > 0 ) ) {
      global.stats.perm <- lapply( stats.across.site.sets.accumulation.fns, function( .fn ) { apply( test.statistics.perm, 2, function( .test.statistics.perm ) { .fn( .test.statistics.perm ) } ) } );
      global.pvals <- lapply( names( global.stats.perm ), function( .global.stat.name ) { calculatePValue.generic( global.stats[[ .global.stat.name ]], global.stats.perm[[ .global.stat.name ]] ) } );
      names( global.pvals ) <- names( global.stats.perm );
    } else {
      global.stats.perm <- NA;
      global.pvals <- NA;
    }
  }

  if( !is.null( site.sets.list ) ) {
    names( pvals ) <- names( site.sets.list );
  }
  # Don't return NA p-values for sites with no variation; instead return 1.
  pvals[ is.na( pvals ) ] <- 1;
  out <- list( p.values=pvals, method = "t.dist")
  if(!use.t.distribution.instead.of.permutation){
    out <- list( p.values=pvals, test.stats=test.statistics, global.pvals=global.pvals, global.stats=global.stats, global.stats.perm=global.stats.perm, include.site=include.site, method=if( use.t.distribution.instead.of.permutation ) { "t.dist" } else if( use.bootstrap ) { "bootstrap" } else { "permutation" } );
  }
  return( out );
} ## tTestPermPval (..)

tTestPermPvalWithRefinement <- function ( vaccine.fasta.file.name, placebo.fasta.file.name, insert.fasta.file.name, weight.matrix.file.name = NULL, include.site = NULL, num.perm = 1000, num.perm.refined = ( num.perm * 10 ), p.value.refine.threshold = .15, be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), ... )
{
  vaccine.tempfile <- NULL;
  if( is.environment( vaccine.fasta.file.name ) ) {
    # Then it will be modified in the first call because of the call to screenOutSites(..).  Let's just write it to a temporary file.
    vaccine.tempfile <- tempfile( fileext = ".fasta" );
    if( be.verbose ) {
      cat( "Saving the given vaccine sequences to a temporary file:", vaccine.tempfile, fill = TRUE );
    }
    writeFastaFile( as.list( vaccine.fasta.file.name ), vaccine.tempfile );
  }
  placebo.tempfile <- NULL;
  if( is.environment( placebo.fasta.file.name ) ) {
    # Then it will be modified in the first call because of the call to screenOutSites(..).  Let's just write it to a temporary file.
    placebo.tempfile <- tempfile( fileext = ".fasta" );
    if( be.verbose ) {
      cat( "Saving the given placebo sequences to a temporary file:", placebo.tempfile, fill = TRUE );
    }
    writeFastaFile( as.list( placebo.fasta.file.name ), placebo.tempfile );
  }
  
  the.result <- tTestPermPval( if( is.null( vaccine.tempfile ) ) { vaccine.fasta.file.name } else { vaccine.tempfile }, if( is.null( placebo.tempfile ) ) { placebo.fasta.file.name } else { placebo.tempfile }, insert.fasta.file.name, weight.matrix.file.name=weight.matrix.file.name, include.site=include.site, num.perm=num.perm, be.verbose = be.verbose, acceptable.chars=acceptable.chars);
   include.site.refined <- the.result$include.site;
   names( the.result$p.values ) <-  which( include.site.refined );
   include.site.refined[ as.numeric( names( the.result$p.values[ the.result$p.values > p.value.refine.threshold ] ) ) ] <- FALSE;
   if( ( num.perm.refined > 0 ) && any( include.site.refined ) ) {
     if( be.verbose ) {
       cat( "Refining ", sum( include.site.refined ), " positions: ", paste( which( include.site.refined ), collapse = ", " ), sep = "", fill = TRUE );
     }
     the.result.refined <- tTestPermPval( if( is.null( vaccine.tempfile ) ) { vaccine.fasta.file.name } else { vaccine.tempfile }, if( is.null( placebo.tempfile ) ) { placebo.fasta.file.name } else { placebo.tempfile }, insert.fasta.file.name, weight.matrix.file.name=weight.matrix.file.name, num.perm=num.perm.refined, include.site=include.site.refined, be.verbose = be.verbose, acceptable.chars=acceptable.chars);
     names( the.result.refined$p.values ) <-  which( include.site.refined );
     the.result$p.values[ as.character( names( the.result.refined$p.values ) ) ] <- the.result.refined$p.values;
   }
  if( !is.null( vaccine.tempfile ) ) {
    unlink( vaccine.tempfile );
  }
  if( !is.null( placebo.tempfile ) ) {
    unlink( placebo.tempfile );
  }
  return( the.result );
} # tTestPermPvalWithRefinement (..)

tTestPermPvalWithRefinement.verifySizeControl <- function ( vaccine.fasta.file.name, placebo.fasta.file.name, insert.fasta.file.name, weight.matrix.file.name = NULL, include.site = NULL, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, num.perm.verifySizeControl = 1000, be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), ... )
{
  # The "file.name" arguments can also take environments.  So our trick is to create permuted environments..
  if( is.environment( insert.fasta.file.name ) ) {
    insert.env <- insert.fasta.file.name;
  } else {
    insert.env <- readFastaFile( insert.fasta.file.name, be.verbose=be.verbose, .acceptable.characters=acceptable.chars );
  }
  if( is.environment( vaccine.fasta.file.name ) ) {
    vaccine.env <- vaccine.fasta.file.name;
  } else {
    vaccine.env <- readFastaFile( vaccine.fasta.file.name, be.verbose=be.verbose, .acceptable.characters=acceptable.chars );
  }
  if( is.environment( placebo.fasta.file.name ) ) {
    placebo.env <- placebo.fasta.file.name;
  } else {
    placebo.env <- readFastaFile( placebo.fasta.file.name, be.verbose=be.verbose, .acceptable.characters=acceptable.chars );
  }

  all.sequences.list <- c( as.list( placebo.env ), as.list( vaccine.env ) );
  
  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.env ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.env );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    placebo.sequence.sets.list = lapply( ls( env=placebo.env ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.env );
  }

  all.sequence.sets.list <- c( vaccine.sequence.sets.list, placebo.sequence.sets.list );
  
  vaccine.num.sequence.sets <- length( vaccine.sequence.sets.list );
  placebo.num.sequence.sets <- length( placebo.sequence.sets.list );
  num.sequence.sets <- vaccine.num.sequence.sets + placebo.num.sequence.sets;
  stopifnot( num.sequence.sets == length( all.sequence.sets.list ) );

  verifySizeControl.p.values <- NULL;
  verifySizeControl.q.values <- NULL;
  for( .perm.i in 1:num.perm.verifySizeControl ) {
    if( be.verbose ) {
      cat( paste( rep( "=", 80 ), collapse = "" ), "\n" );
      cat( paste( "[tTestPermPvalWithRefinement.verifySizeControl(..)] Calling tTestPermPvalWithRefinement(..) for permutation", .perm.i, "of", num.perm.verifySizeControl, "..\n" ) );
    }
    .permutation <- sample( 1:num.sequence.sets, replace = FALSE );

    .perm.vaccine.sequence.sets.list <- all.sequence.sets.list[ .permutation[ 1:vaccine.num.sequence.sets ] ];
    .perm.placebo.sequence.sets.list <- all.sequence.sets.list[ .permutation[ ( vaccine.num.sequence.sets + 1 ):num.sequence.sets ] ];

    .perm.vaccine.env <- new.env( hash = TRUE );
    lapply( .perm.vaccine.sequence.sets.list, function( .sequence.set ) { sapply( .sequence.set, function( .sequence.name ) { assign( .sequence.name, all.sequences.list[[ .sequence.name ]], env=.perm.vaccine.env ) } ) } );
    .perm.placebo.env <- new.env( hash = TRUE );
    lapply( .perm.placebo.sequence.sets.list, function( .sequence.set ) { sapply( .sequence.set, function( .sequence.name ) { assign( .sequence.name, all.sequences.list[[ .sequence.name ]], env=.perm.placebo.env ) } ) } );
    
    .result <- tTestPermPvalWithRefinement( .perm.vaccine.env, .perm.placebo.env, insert.fasta.file.name, weight.matrix.file.name=weight.matrix.file.name, include.site=include.site, vaccine.sequence.sets.list=.perm.vaccine.sequence.sets.list, placebo.sequence.sets.list=.perm.placebo.sequence.sets.list, be.verbose = be.verbose, ... );
    if( .perm.i == 1 ) {
      verifySizeControl.p.values <- matrix( nrow = num.perm.verifySizeControl, ncol = length( .result$p.values ) );
      verifySizeControl.q.values <- matrix( nrow = num.perm.verifySizeControl, ncol = length( .result$p.values ) );
    }
    verifySizeControl.p.values[ .perm.i, ] <- .result$p.values;
    verifySizeControl.q.values[ .perm.i, ] <- p.adjust( .result$p.values, method = "fdr" );
  } # End foreach .perm.i
  return( list( p.values = verifySizeControl.p.values, q.values = verifySizeControl.q.values ) );
} # tTestPermPvalWithRefinement.verifySizeControl (..)

######################### Non-insert-dependent ###############
### 'klStat' returns values of the Kullback-Leibler test statistic per aa position.
### NOTE: this is NOT symmetric!  It calculates KL( V || P ), which is different from KL( P || V ).  If you want a symmeterized kl divergence, you can use sklStat(..) instead.  See also hellingerStat(..).
#### TODO: Adapt to take acceptable.chars, to handle mising data.
klStat <- function ( vaccine.seq.envir, placebo.seq.envir ) {
	vaccine.results <- aaProportionsPerPosition(vaccine.seq.envir)
	placebo.results <- aaProportionsPerPosition(placebo.seq.envir)

	if (vaccine.results$seq.length != placebo.results$seq.length){ stop("Vaccine and placebo sequences are of unequal length.") }
	
	vaccine.props <- vaccine.results$aa.proportions
	placebo.props <- placebo.results$aa.proportions

	num.vaccine.seqs <- vaccine.results$num.seqs
	num.placebo.seqs <- placebo.results$num.seqs

	seq.length <- vaccine.results$seq.length

	test.statistics <- numeric(seq.length)
	discrepancy <- function(a, b){ a * log(a/b) }
	for (i in 1:seq.length){
		test.statistics[i] <- sum( ifelse(vaccine.props[[i]]*placebo.props[[i]]>0, discrepancy(vaccine.props[[i]], placebo.props[[i]]), discrepancy(vaccine.props[[i]] + 1/num.vaccine.seqs, placebo.props[[i]] + 1/num.placebo.seqs)) )
	}
	return( test.statistics )
} ## klStat (..)

sklStat <- function ( vaccine.seq.envir, placebo.seq.envir ) {
  return( klStat( vaccine.seq.envir, placebo.seq.envir ) + klStat( placebo.seq.envir, vaccine.seq.envir ) );
} ## sklStat (..)

### return the hellinger distances at each position across the two groups.
hellingerStat <- function ( vaccine.seq.envir, placebo.seq.envir ) {
  vaccine.proportionsPerPosition <- aaProportionsPerPosition( vaccine.seq.envir );
  placebo.proportionsPerPosition <- aaProportionsPerPosition( placebo.seq.envir );

  if( vaccine.proportionsPerPosition$seq.length != placebo.proportionsPerPosition$seq.length ){
    stop( "Vaccine and placebo sequences are of unequal length." );
  }
	
  vaccine.props <- vaccine.proportionsPerPosition$aa.proportions;
  placebo.props <- placebo.proportionsPerPosition$aa.proportions;

  num.vaccine.seqs <- vaccine.proportionsPerPosition$num.seqs;
  num.placebo.seqs <- placebo.proportionsPerPosition$num.seqs;

  seq.length <- vaccine.proportionsPerPosition$seq.length;
  
  test.statistics <- numeric( seq.length );
  for( i in 1:seq.length ) {
    test.statistics[ i ] <- sqrt( 1 - sum( sqrt( vaccine.props[[i]]*placebo.props[[i]] ) ) );
  }
  return( test.statistics )
} # hellingerStat (..)

### 'klPermTestStatistics' returns a matrix of Kullback-Leibler (or SKL, or Hellinger) test statistics per aa position from 'num.perm' permuted data sets.
## Set use.skl.test to TRUE to use "sklStat(..)" instead of "klStat".
## Set use.hellinger.test to TRUE to use "hellingerStat(..)" instead of "klStat".
klPermTestStatistics <- function ( vaccine.seq.env, placebo.seq.env, num.perm, use.skl.test = FALSE, use.hellinger.test = FALSE, be.verbose = FALSE )
{
  if( use.skl.test & use.hellinger.test ) {
    stop( "You may specify that klPermTestStatistics should use.skl.test or that it should use.hellinger.test, but not both." );
  }
	all.seq.list <- as.list(c(do.call("c", as.list(vaccine.seq.env)), do.call("c", as.list(placebo.seq.env))))
	vaccine.num.seqs <- countFastaSeqs(vaccine.seq.env)
	placebo.num.seqs <- countFastaSeqs(placebo.seq.env)
	num.seqs <- vaccine.num.seqs + placebo.num.seqs
	seq.length <- nchar(get(ls(name=vaccine.seq.env)[1], env=vaccine.seq.env))

	test.statistics.perm <- matrix(0, nrow=seq.length, ncol=num.perm)
	for (B in 1:num.perm){
          if( be.verbose && ( ( B %% floor( num.perm / 100 ) ) == 0 ) ) {
            if( use.hellinger.test ) {
              cat( paste( "Calling hellingerStat(..) for permutation", B, "of", num.perm, "..\n" ) );
            } else if( use.skl.test ) {
              cat( paste( "Calling sklStat(..) for permutation", B, "of", num.perm, "..\n" ) );
            } else {
              cat( paste( "Calling klStat(..) for permutation", B, "of", num.perm, "..\n" ) );
            }
          }
		permutation <- sample(1:num.seqs)

		vaccine.perm.env <- new.env(hash=TRUE)
		for (j in permutation[1:vaccine.num.seqs]){ assign(names(all.seq.list)[j], all.seq.list[[j]], env=vaccine.perm.env) }

		placebo.perm.env <- new.env(hash=TRUE)
		for (j in permutation[(vaccine.num.seqs+1):num.seqs]){ assign(names(all.seq.list)[j], all.seq.list[[j]], env=placebo.perm.env) }
                if( use.hellinger.test ) {
                  test.statistics.perm[ , B ] <- hellingerStat( vaccine.perm.env, placebo.perm.env );
                } else if( use.skl.test ) {
                  test.statistics.perm[ , B ] <- sklStat( vaccine.perm.env, placebo.perm.env );
                } else {
                  test.statistics.perm[ , B ] <- klStat( vaccine.perm.env, placebo.perm.env );
                }
	}
	return( test.statistics.perm )
} # klPermTestStatistics (..)

sklPermTestStatistics <- function ( ... )
{
  return( klPermTestStatistics( ..., use.skl.test = TRUE ) );
} # sklPermTestStatistics (..)

hellingerPermTestStatistics <- function ( ... )
{
  return( klPermTestStatistics( ..., use.hellinger.test = TRUE ) );
} # hellingerPermTestStatistics (..)

##### NOTE: the kl test is NOT the skl test.  That is, it matters which group is vax, which is pla.  However, if you specify use.skl.test = TRUE, then it will be the skl test instead.  If you instead specify use.hellinger.test = TRUE, then that will be used.  Don't try to use both simultaneously, though.
### Note that sklTestPermPval(..) is a convenience method that delegates to this after setting use.skl.test to TRUE.
### Note that hellingerTestPermPval(..) is a convenience method that delegates to this after setting use.hellinger.test to TRUE.
### 'klTestPermPval' returns a list with the following components:
### 'p.values'               - numeric vector of unadjusted permutation-based Kullback-Leibler (or SKL, or Hellinger) test p-values per aa position
### 'significance'           - character vector with strings "reject" and "do not reject" to determine significance
### 'significance.threshold' - fine tuned p.cut to match the estimated and nominal numbers of FP rejections
klTestPermPval <- function ( vaccine.fasta.file.name, placebo.fasta.file.name, include.site=NULL, num.perm=100, use.skl.test = FALSE, use.hellinger.test = FALSE, determine.signif=FALSE, num.perm.signif=NULL, num.FP.rejections=NULL, be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{
  if( use.skl.test & use.hellinger.test ) {
    stop( "You may specify that klTestPermPval should use.skl.test or that it should use.hellinger.test, but not both." );
  }

  if( is.environment( vaccine.fasta.file.name ) ) {
    vaccine.env <- vaccine.fasta.file.name;
  } else {
    vaccine.env <- readFastaFile(vaccine.fasta.file.name, .acceptable.characters=acceptable.chars);
  }

  if( is.environment( placebo.fasta.file.name ) ) {
    placebo.env <- placebo.fasta.file.name;
  } else {
    placebo.env <- readFastaFile(placebo.fasta.file.name, .acceptable.characters=acceptable.chars);
  }

	if (!is.null(include.site)){
		vaccine.env <- screenOutSites(vaccine.env, include.site)
		placebo.env <- screenOutSites(placebo.env, include.site)
	}

	all.seq.list <- as.list(c(do.call("c", as.list(vaccine.env)), do.call("c", as.list(placebo.env))))
	vaccine.num.seqs <- countFastaSeqs(vaccine.env)
	placebo.num.seqs <- countFastaSeqs(placebo.env)
	num.seqs <- vaccine.num.seqs + placebo.num.seqs

        if( use.hellinger.test ) {
          if( be.verbose ) {
            cat( "Calling hellingerStat(..) for main test stats..\n" );
          }
          test.statistics <- hellingerStat( vaccine.env, placebo.env );
          test.statistics.perm <- hellingerPermTestStatistics( vaccine.env, placebo.env, num.perm, be.verbose=be.verbose );
        } else if( use.skl.test ) {
          if( be.verbose ) {
            cat( "Calling sklStat(..) for main test stats..\n" );
          }
          test.statistics <- sklStat( vaccine.env, placebo.env );
          test.statistics.perm <- sklPermTestStatistics( vaccine.env, placebo.env, num.perm, be.verbose=be.verbose );
        } else {
          if( be.verbose ) {
            cat( "Calling klStat(..) for main test stats..\n" );
          }
          test.statistics <- klStat( vaccine.env, placebo.env );
          test.statistics.perm <- klPermTestStatistics( vaccine.env, placebo.env, num.perm, be.verbose=be.verbose );
        }
        ##### NOTE: The skl test stat is NOT symmetric about 0 under the null.  Nor is hellinger.  I think that the kl test should be symmetric about 0 but I've not verified it [I believe its expectation is 0, so at least it's centered at 0, but I worry that unequal sample sizes might make it non-symmetric even under the null].  Thus we do not use this:
	#pvals <- apply(abs(test.statistics.perm) >= abs(test.statistics), 1, mean)
        ##### NOTE: instead we use this:
        calculatePValue <- function ( .test.statistic, .null.test.statistics, add.one.for.observed.test.statistic = FALSE ) { min( 1, 2 * ( min( sum( .test.statistic >= .null.test.statistics, na.rm=T ), sum( .test.statistic <= .null.test.statistics, na.rm=T ) ) + as.numeric( add.one.for.observed.test.statistic ) ) / (  sum( !is.na( .null.test.statistics ) ) + as.numeric( add.one.for.observed.test.statistic ) ), na.rm=T ) }; 
        pvals <- sapply( 1:nrow( test.statistics.perm ), function( pos.i ) { calculatePValue( test.statistics[ pos.i ], test.statistics.perm[ pos.i, ] ) } );

	if (determine.signif){
		### Step (1) in GWJ Section 3.3 ###
		new.permutation <- sample(1:num.seqs)

		vaccine.new.env <- new.env(hash=TRUE)
		for (j in new.permutation[1:vaccine.num.seqs]){ assign(names(all.seq.list)[j], all.seq.list[[j]], env=vaccine.new.env) }

		placebo.new.env <- new.env(hash=TRUE)
		for (j in new.permutation[(vaccine.num.seqs+1):num.seqs]){ assign(names(all.seq.list)[j], all.seq.list[[j]], env=placebo.new.env) }

		### Step (2) in GWJ Section 3.3 ###
                if( use.hellinger.test ) {
                  test.statistics.signif <- hellingerStat( vaccine.new.env, placebo.new.env );
                  test.statistics.perm.signif <- hellingerPermTestStatistics( vaccine.new.env, placebo.new.env, num.perm.signif );
                } else if( use.skl.test ) {
                  test.statistics.signif <- sklStat( vaccine.new.env, placebo.new.env );
                  test.statistics.perm.signif <- sklPermTestStatistics( vaccine.new.env, placebo.new.env, num.perm.signif );
                } else {
                  test.statistics.signif <- klStat( vaccine.new.env, placebo.new.env );
                  test.statistics.perm.signif <- klPermTestStatistics( vaccine.new.env, placebo.new.env, num.perm.signif );
                }
                ### SEE NOTE ABOVE about this calculation.
                pvals.signif <- sapply( 1:nrow( test.statistics.perm.signif ), function( pos.i ) { calculatePValue( test.statistics.signif[ pos.i ], test.statistics.perm.signif[ pos.i, ] ) } );

		### Steps (3) and (4) in GWJ Section 3.3 ###
		f <- function(p.cut, pval.vect, num.FP.rejections){ sum(pval.vect < p.cut) - num.FP.rejections }
		p.cut.tuned <- uniroot(f, 0:1, pval.vect=pvals.signif, num.FP.rejections=num.FP.rejections)$root

		### Step (5) in GWJ Section 3.3 ###
		rejection.results <- ifelse(pvals < p.cut.tuned, "reject", "do not reject")
	}
	out <- list(p.values=pvals)
	if (determine.signif){
		out$significance <- rejection.results
		out$significance.threshold <- p.cut.tuned
	}
	return( out )
} ## klTestPermPval(..)

sklTestPermPval <- function ( ... ) {
  return( klTestPermPval( ..., use.skl.test = TRUE ) );
} ## sklTestPermPval (..)

hellingerTestPermPval <- function ( ... ) {
  return( klTestPermPval( ..., use.hellinger.test = TRUE ) );
} ## hellingerTestPermPval (..)


# The following 6 functions were appended to this file from Mullins454STEP_sieve_454_methods_mba.R on 7/12/2013
#   1.  tStatbyPos.454 
#   2.	tPermTestStatisticsbyPos.454 
#   3.	tTestPermPval.454
#   4.	tTestPermPvalWithRefinement.454
#   5.	tTestPermPvalWithRefinement.verifySizeControl.454
#   6.  chooseSitesToScreen.minimumVariability.454

#############################################tStatbyPos.454#####################################################
# If min.variant.freq.threshold > 0, first set any frequencies <= min.variant.freq.threshold to 0, and then renormalize, before calculating weights.
tStatbyPos.454 <- function( vaccine.aa.freq.list.pos, placebo.aa.freq.list.pos, insert.char.pos, weight.matrix = NULL, return.t.test.result = FALSE, acceptable.chars = if( is.null( weight.matrix )){ c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) } else{ setdiff( rownames( weight.matrix ), 'X' ) }, min.variant.freq.threshold = 0, normalize.frequencies = TRUE, return.list = FALSE )
{
  if( is.na( insert.char.pos ) ) {
    return( NA );
  }
  # make sure the insert char is an allowable character.. 
  stopifnot( insert.char.pos %in% acceptable.chars );
  
  vaccine.num.sequences <- length( vaccine.aa.freq.list.pos );
  vaccine.weights <- rep( NA, vaccine.num.sequences );
  names( vaccine.weights ) <- names( vaccine.aa.freq.list.pos ); # assumes rownames are sequence names
  for( .subject.k in 1:length( vaccine.weights ) ) {
    .vaccine.aa.freqs.for.subject.at.pos <- vaccine.aa.freq.list.pos[[ .subject.k ]];
    # Replace values <= min.variant.freq.threshold with 0.
    .vaccine.aa.freqs.for.subject.at.pos[ .vaccine.aa.freqs.for.subject.at.pos <= min.variant.freq.threshold ] <- 0;
    if( normalize.frequencies ) {
      .vaccine.aa.freqs.for.subject.at.pos <- .vaccine.aa.freqs.for.subject.at.pos / sum( .vaccine.aa.freqs.for.subject.at.pos, na.rm = TRUE );
    }
    if( all( is.na( .vaccine.aa.freqs.for.subject.at.pos ) ) ) {  # NOTE: assumes subjects with entirely-missing data are represented by a row of NAs
      .weight <- NA;
    } else if( is.null( weight.matrix ) ) {
      # Then use a 0-1 weight matrix (so the .weight is just the freq of the insert char).
      .weight <- .vaccine.aa.freqs.for.subject.at.pos[ insert.char.pos ];
    } else {
      stopifnot( !any( is.na( .vaccine.aa.freqs.for.subject.at.pos ) ) ); # See above; if ANY are NA, then ALL should be.
      # take the expectation of weight across AA distribution
      .weight <- .vaccine.aa.freqs.for.subject.at.pos %*% weight.matrix[ insert.char.pos, ];
    }
    vaccine.weights[ .subject.k ] <- .weight;
  } # End foreach vaccine seq.k

  placebo.num.sequences <- length( placebo.aa.freq.list.pos );
  placebo.weights <- rep( NA, placebo.num.sequences );
  names( placebo.weights ) <- names( placebo.aa.freq.list.pos ); # assumes rownames are sequence names
  for( .subject.k in 1:length( placebo.weights ) ) {
    .placebo.aa.freqs.for.subject.at.pos <- placebo.aa.freq.list.pos[[ .subject.k ]];
    # Replace values <= min.variant.freq.threshold with 0.
    .placebo.aa.freqs.for.subject.at.pos[ .placebo.aa.freqs.for.subject.at.pos <= min.variant.freq.threshold ] <- 0;
    if( normalize.frequencies ) {
      .placebo.aa.freqs.for.subject.at.pos <- .placebo.aa.freqs.for.subject.at.pos / sum( .placebo.aa.freqs.for.subject.at.pos, na.rm = TRUE );
    }
    if( all( is.na( .placebo.aa.freqs.for.subject.at.pos ) ) ) {  # NOTE: assumes subjects with entirely-missing data are represented by a row of NAs
      .weight <- NA;
    } else if( is.null( weight.matrix ) ) {
      # Then use a 0-1 weight matrix (so the .weight is just the freq of the insert char).
      .weight <- .placebo.aa.freqs.for.subject.at.pos[ insert.char.pos ];
    } else {
      stopifnot( !any( is.na( .placebo.aa.freqs.for.subject.at.pos ) ) ); # See above; if ANY are NA, then ALL should be.
      # take the expectation of weight across AA distribution
      .weight <- .placebo.aa.freqs.for.subject.at.pos %*% weight.matrix[ insert.char.pos, ];
    }
    placebo.weights[ .subject.k ] <- .weight;
  } # End foreach placebo seq.k

  test.statistic <- NA;
  
  num.vaccine.weights <- sum( !is.na( vaccine.weights ) );
  num.placebo.weights <- sum( !is.na( placebo.weights ) );
  num.weights <- num.vaccine.weights + num.placebo.weights;
  
  if( return.t.test.result ) {
    if( all( vaccine.weights == vaccine.weights[ 1 ] ) && all( placebo.weights == vaccine.weights[ 1 ] ) ) {
      # Then there's no variation.
      test.statistic <- NA;
    } else {
      t.test.result <- t.test( vaccine.weights, placebo.weights );
      test.statistic <- list( t.test.result = t.test.result );
    }
  } else {
    se <- sqrt( (num.vaccine.weights - 1)*var(vaccine.weights, na.rm=T)/(num.weights - 2) + (num.placebo.weights - 1)*var(placebo.weights, na.rm=T)/(num.weights - 2) );
    test.statistic <- ( ( mean( vaccine.weights, na.rm=T ) - mean( placebo.weights, na.rm=T ) ) / se );
  }
  if( return.list ) {
    return( test.statistic = test.statistic, vaccine.weights = vaccine.weights, placebo.weights = placebo.weights, se = se );
  } else {
    return( test.statistic );
  }
} ## tStatbyPos.454 (..)

############################################tPermTestStatisticsbyPos.454####################################################
tPermTestStatisticsbyPos.454 <- function( vaccine.aa.freq.list.pos, placebo.aa.freq.list.pos, insert.char.pos, weight.matrix = NULL, return.t.test.result = FALSE, acceptable.chars = if( is.null( weight.matrix )){ c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) } else{ setdiff( rownames( weight.matrix ), 'X' ) }, num.perm = 100, use.bootstrap = FALSE, be.verbose = TRUE, min.variant.freq.threshold = 0, normalize.frequencies = TRUE, return.list = FALSE )
{
  combined.aa.freq.list.pos <- c( vaccine.aa.freq.list.pos, placebo.aa.freq.list.pos );
  
  vaccine.num.sequences <- length( vaccine.aa.freq.list.pos );
  placebo.num.sequences <- length( placebo.aa.freq.list.pos );
  num.sequences <- vaccine.num.sequences + placebo.num.sequences;
  
  stopifnot( num.sequences == length( combined.aa.freq.list.pos ) );
  test.statistics.perm <- rep( 0, times = num.perm );
  if( return.list ) {
    test.results.list.perm <- list();
  }
  for( perm.i in 1:num.perm ) {
    if( be.verbose && ( ( num.perm <= 100 ) || ( ( perm.i %% floor( num.perm / 100 ) ) == 0 ) ) ) {
      cat( paste( "Calling tStatbyPos.454(..) for permutation", perm.i, "of", num.perm, "..\n" ) );
    }
    .permutation <- sample( 1:num.sequences, replace = use.bootstrap );
    .perm.vaccine.aa.freq.list.pos <- combined.aa.freq.list.pos[ .permutation[ 1:vaccine.num.sequences ] ];
    .perm.placebo.aa.freq.list.pos <- combined.aa.freq.list.pos[ .permutation[ ( vaccine.num.sequences + 1 ):num.sequences ] ];
    .test.statistics <- tStatbyPos.454( vaccine.aa.freq.list.pos=.perm.vaccine.aa.freq.list.pos, placebo.aa.freq.list.pos=.perm.placebo.aa.freq.list.pos, insert.char.pos=insert.char.pos, weight.matrix=weight.matrix, return.t.test.result=return.t.test.result, acceptable.chars=acceptable.chars, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list );
    if( return.list ) {
      .test.results.list <- .test.statistics;
      test.results.list.perm[[ perm.i ]] <- .test.results.list;
      .test.statistics <- .test.statistics$test.statistics;
    }
    test.statistics.perm[ perm.i ] <- .test.statistics;
  } # End foreach permutation perm.i
  if( return.list ) {
    return( list( test.statistics.perm = test.statistics.perm, test.results.list.perm = test.results.list.perm ) );
  } else {
    return( test.statistics.perm );
  }
} # tPermTestStatisticsbyPos.454 (..)

########################################################tTestPermPval.454##################################################################
tTestPermPval.454 <- function ( vaccine.aa.freq.list, placebo.aa.freq.list, insert.fasta.file.name, weight.matrix.file.name = NULL, weight.matrix.is.probability.matrix = FALSE, include.site = NULL, num.perm = 100, use.t.distribution.instead.of.permutation=( num.perm == 0 ), use.bootstrap =FALSE, be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), min.variant.freq.threshold = 0, normalize.frequencies = TRUE, return.list = FALSE, ... )
{
  if( is.null( weight.matrix.file.name ) ) {
    # If weight.matrix is NULL, weights will be 0 for insert-match and 1 for insert-mismatch.
    weight.matrix <- NULL;
    if( be.verbose ) {
      cat( "Using mismatch indicators instead of weights.", fill = T );
    }
  } else {
    if( weight.matrix.is.probability.matrix ) {
      # In probability form
      if( is.matrix( weight.matrix.file.name ) ) {
        weight.matrix <- 1 - weight.matrix.file.name;   ## convert to a dissimilarity matrix
      } else {
        weight.matrix <- 1 - readBetty(weight.matrix.file.name);   ## convert to a dissimilarity matrix
      }
    } else {
      # In log-likelihood-ratio (etc) form
      if( is.matrix( weight.matrix.file.name ) ) {
        weight.matrix <- 0 - weight.matrix.file.name;    ## convert to a dissimilarity matrix
      } else {
        weight.matrix <- 0 - readPAM(weight.matrix.file.name);    ## convert to a dissimilarity matrix
      }
    }
    if( be.verbose ) {
      cat( "Using dissimilarity matrix converted from the", ifelse( weight.matrix.is.probability.matrix, "probability-form substitution matrix", "substitution matrix" ), "in file:", weight.matrix.file.name, fill = T );
    }
  } # End if is.null( weight.matrix.file.name ) .. else ..
  
  if( is.environment( insert.fasta.file.name ) ) {
    insert.env <- insert.fasta.file.name;
  } else {
    insert.env <- readFastaFile( insert.fasta.file.name, be.verbose=be.verbose, .acceptable.characters = acceptable.chars );
  }
  insert.char.vector <- strsplit(as.list(insert.env)[[1]],"")[[1]];
  
  ## TODO: ? .. For now, if include.site is _any_ character string, then use the minimumVariabilty rules.
  if( is.character( include.site ) ) {
    if( be.verbose ) {
      cat( "Calling chooseSitesToScreen.minimumVariability.454(..)", fill = T );
    }
    include.site <- chooseSitesToScreen.minimumVariability.454( vaccine.aa.freq.list=vaccine.aa.freq.list, placebo.aa.freq.list=placebo.aa.freq.list, insert.seq.envir=insert.env, min.variant.freq.threshold=min.variant.freq.threshold, ... );
  }
  include.sites <- which( include.site == TRUE );
  
  ## NOTE we're excluding the site.set and multiple-seqs-per-subject list 
  
  ### NOTE that these are two-sided tests!
  if( use.t.distribution.instead.of.permutation ) {
    test.statistics <- sapply( include.sites, function( .pos ){ tStatbyPos.454( vaccine.aa.freq.list.pos = vaccine.aa.freq.list[[ .pos ]], placebo.aa.freq.list.pos = placebo.aa.freq.list[[ .pos ]], insert.char.pos = insert.char.vector[ .pos ], weight.matrix = weight.matrix, return.t.test.result = TRUE, acceptable.chars = acceptable.chars, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list ); } );
    if( return.list ) {
      test.results.list <- test.statistics;
      test.statistics <- test.statistics$test.statistics;
    }
    #pvals <- sapply( test.statistics, function( .t.test.result ){ if( is.na( .t.test.result ) ) { return( NA ); } else { return( .t.test.result$p.value ) } } );
    pvals <- sapply( test.statistics, function( .t.test.result ) return (ifelse(is.null(as.list(.t.test.result)$p.value), .t.test.result, .t.test.result$p.value ))  )
    } else {
    if( be.verbose ) {
      cat( "Calling tStatbyPos.454(..) for main test stats..\n" );
    }
    test.statistics <- sapply( include.sites, function( .pos ) { tStatbyPos.454( vaccine.aa.freq.list.pos = vaccine.aa.freq.list[[ .pos ]], placebo.aa.freq.list.pos = placebo.aa.freq.list[[ .pos ]], insert.char.pos = insert.char.vector[ .pos ], weight.matrix = weight.matrix, return.t.test.result = FALSE, acceptable.chars = acceptable.chars, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list ); } );
    if( return.list ) {
      test.results.list <- test.statistics;
      test.statistics <- test.statistics$test.statistics;
    }
    if( be.verbose ) {
      cat( "Calling tPermTestStatisticsbyPos.454(..) for", ifelse( use.bootstrap, "bootstraps", "permutations" ), "test stats..\n" );
    }
    test.statistics.perm <- t( as.matrix( sapply( include.sites, function( .pos ){ tPermTestStatisticsbyPos.454( vaccine.aa.freq.list.pos = vaccine.aa.freq.list[[ .pos ]], placebo.aa.freq.list.pos = placebo.aa.freq.list[[ .pos ]], insert.char.pos = insert.char.vector[ .pos ], weight.matrix = weight.matrix, num.perm = num.perm, use.bootstrap = use.bootstrap, be.verbose = be.verbose, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list ) } ) ) );
    if( return.list ) {
      test.results.list <- test.statistics.perm;
      test.statistics.perm <- test.results.list$test.statistics.perm;
    }
    ##### NOTE: The tTest stat should be symmetric about 0 under the null but I've not verified it.  Thus we do not use this:
    #pvals <- apply( abs(test.statistics.perm) >= abs(test.statistics), 1, mean );
    ##### NOTE: instead we use this:
    calculatePValue <- function ( .pos.i, add.one.for.observed.test.statistic = FALSE ) { if( is.na( test.statistics[ .pos.i ] ) ) { return( NA ); } else { return( min( 1, 2 * ( min( sum( test.statistics[ .pos.i ] >= test.statistics.perm[ .pos.i,  ], na.rm=T ), sum( test.statistics[ .pos.i ] <= test.statistics.perm[ .pos.i,  ], na.rm=T ) ) + as.numeric( add.one.for.observed.test.statistic ) ) / (  sum( !is.na( test.statistics.perm[ .pos.i,  ] ) ) + as.numeric( add.one.for.observed.test.statistic ) ), na.rm=T ) ); } };
    pvals <- sapply( 1:nrow( test.statistics.perm ), calculatePValue );
    ## TODO: REMOVE
    if( any( is.na( pvals ) ) ) {
      print( "WEIRD: GOT NA pval for eg" )
      na.pval.pos <- which( is.na( pvals ) )[ 1 ];
      print( na.pval.pos );
      print( test.statistics[ na.pval.pos ] );
      print( test.statistics.perm[ na.pval.pos, ] );
    }
  }
  if( !is.null( include.sites ) ) {
    names( pvals ) <- include.sites;
  }
  out <- list( p.values=pvals, test.stats=test.statistics, include.site = include.site, method=if( use.t.distribution.instead.of.permutation ) { "t.dist" } else if( use.bootstrap ) { "bootstrap" } else { "permutation" } );
  ## TODO: REMOVE.  Including test.statistics.perm for debugging.
  #out[[ "test.statistics.perm" ]] <- test.statistics.perm;
  if( return.list ) {
    out[[ "test.results.list" ]] <- test.results.list;
  }
  return( out );
} ## tTestPermPval.454 (..)

########################################tTestPermPvalWithRefinement.454#################################
tTestPermPvalWithRefinement.454 <- function ( vaccine.aa.freq.list, placebo.aa.freq.list, insert.fasta.file.name, weight.matrix.file.name = NULL, include.site = NULL, num.perm = 1000, num.perm.refined = ( num.perm * 10 ), p.value.refine.threshold = .15, be.verbose = TRUE, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), min.variant.freq.threshold = 0, normalize.frequencies = TRUE, return.list = FALSE, ... )
{
  the.result <- tTestPermPval.454( vaccine.aa.freq.list=vaccine.aa.freq.list, placebo.aa.freq.list=placebo.aa.freq.list, insert.fasta.file.name=insert.fasta.file.name, weight.matrix.file.name=weight.matrix.file.name, include.site=include.site, num.perm=num.perm, be.verbose = be.verbose, acceptable.chars=acceptable.chars, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list, ... );
  include.site.refined <- the.result$include.site;
  names( the.result$p.values ) <-  which( include.site.refined );
  include.site.refined[ as.numeric( names( the.result$p.values[ the.result$p.values > p.value.refine.threshold ] ) ) ] <- FALSE;
  if( ( num.perm.refined > 0 ) && any( include.site.refined ) ) {
    if( be.verbose ) {
      cat( "Refining ", sum( include.site.refined ), " positions: ", paste( which( include.site.refined ), collapse = ", " ), sep = "", fill = TRUE );
    }
    the.result.refined <- tTestPermPval.454( vaccine.aa.freq.list=vaccine.aa.freq.list, placebo.aa.freq.list=placebo.aa.freq.list, insert.fasta.file.name=insert.fasta.file.name, weight.matrix.file.name=weight.matrix.file.name, include.site=include.site.refined, num.perm=num.perm.refined, be.verbose = be.verbose, acceptable.chars=acceptable.chars, min.variant.freq.threshold = min.variant.freq.threshold, normalize.frequencies = normalize.frequencies, return.list = return.list, ... );
    names( the.result.refined$p.values ) <-  which( include.site.refined );
    the.result$p.values[ as.character( names( the.result.refined$p.values ) ) ] <- the.result.refined$p.values;
  }
  return( the.result );
} # tTestPermPvalWithRefinement.454 (..)

################################tTestPermPvalWithRefinement.verifySizeControl.454############################
### TODO: UPDATE/FIX.
# tTestPermPvalWithRefinement.verifySizeControl.454 <- function ( vaccine.aa.freq.list, placebo.aa.freq.list, insert.fasta.file.name, weight.matrix.file.name = NULL, include.site = NULL, num.perm.verifySizeControl = 1000, be.verbose = TRUE, min.variant.freq.threshold, use.aa.freq.sum, ... )
# {
#   # The "file.name" arguments can also take environments.  So our trick is to create permuted environments..
#   if( is.environment( insert.fasta.file.name ) ) {
#     insert.env <- insert.fasta.file.name;
#   } else {
#     insert.env <- readFastaFile( insert.fasta.file.name, be.verbose=be.verbose );
#   }
#   
#   # get the list of all the names:
#   all.sequences.list <- c( as.list( dimnames( vaccine.aa.freq.list[[1]] )[[1]] ), as.list( dimnames( placebo.aa.freq.list[[1]] )[[1]] ) );
#   
#   # find the total # of sequences
#   num.vaccine.seqs <- length( dimnames( vaccine.aa.freq.list[[1]] )[[1]] )
#   num.placebo.seqs <- length( dimnames( placebo.aa.freq.list[[1]] )[[1]] )
#   
#   num.sequences <- num.vaccine.seqs + num.placebo.seqs;
#   
#   verifySizeControl.p.values <- NULL;
#   verifySizeControl.q.values <- NULL;
#   for( .perm.i in 1:num.perm.verifySizeControl ) {
#     if( be.verbose ) {
#       cat( paste( rep( "=", 80 ), collapse = "" ), "\n" );
#       cat( paste( "[tTestPermPvalWithRefinement.verifySizeControl(..)] Calling tTestPermPvalWithRefinement(..) for permutation", .perm.i, "of", num.perm.verifySizeControl, "..\n" ) );
#     }
#     .permutation <- sample( 1:num.sequences, replace = FALSE );
#     
#     .perm.vaccine.aa.freq.list <- list();
#     .perm.placebo.aa.freq.list <- list();
#     for( .pos.i in 1:length( vaccine.aa.freq.list ) ){
#       ## NOTE: FIRST WAY, might be slow... faster than second way
#       .combined.aa.array <- rbind( vaccine.aa.freq.list[[ .pos.i ]], placebo.aa.freq.list[[ .pos.i ]] );
#       .perm.vaccine.aa.freq.list[[ .pos.i ]] <- .combined.aa.array[ .permutation[ 1:num.vaccine.seqs ], ];
#       .perm.placebo.aa.freq.list[[ .pos.i ]] <- .combined.aa.array[ .permutation[ ( num.vaccine.seqs + 1 ):num.sequences ], ];
#       ## SECOND WAY:
#       # .vaccine.indices <- .permutation[ 1:num.vaccine.seqs ];
#       # get the indices of sequences from the vaccine group and the placebo group.
#       #.vax.in.vax <- .vaccine.indices[ .vaccine.indices <= num.vaccine.seqs ];
#       #.vax.in.pla <- .vaccine.indices[ .vaccine.indices > num.vaccine.seqs ] - num.vaccine.seqs;
#       #.perm.vaccine.aa.freq.list[[ .pos.i ]] <- rbind( vaccine.aa.freq.list[[ .pos.i ]][ .vax.in.vax, ], placebo.aa.freq.list[[ .pos.i ]][ .vax.in.pla, ] );
#       #.placebo.indices <- .permutation[ (num.vaccine.seqs + 1):num.sequences ];
#       #.pla.in.vax <- .placebo.indices[ .placebo.indices <= num.vaccine.seqs ];
#       #.pla.in.pla <- .placebo.indices[ .placebo.indices > num.vaccine.seqs ] - num.vaccine.seqs;
#       #.perm.vaccine.aa.freq.list[[ .pos.i ]] <- rbind( vaccine.aa.freq.list[[ .pos.i ]][ .pla.in.vax, ], placebo.aa.freq.list[[ .pos.i ]][ .pla.in.pla ] );
#     } # end foreach .pos.i
#     
#     .result <- tTestPermPvalWithRefinement.454( vaccine.aa.freq.list = .perm.vaccine.aa.freq.list, placebo.aa.freq.list = .perm.placebo.aa.freq.list, insert.fasta.file.name = insert.fasta.file.name, weight.matrix.file.name = weight.matrix.file.name, include.site = include.site, num.perm = 1000, p.value.refine.threshold = .15, be.verbose = be.verbose, min.variant.freq.threshold = min.variant.freq.threshold, use.aa.freq.sum = use.aa.freq.sum, ... );
#     if( .perm.i == 1 ){
#       verifySizeControl.p.values <- matrix( nrow = num.perm.verifySizeControl, ncol = length( .result$p.values ) );
#       verifySizeControl.q.values <- matrix( nrow = num.perm.verifySizeControl, ncol = length( .result$p.values ) );      
#     }
#     verifySizeControl.p.values[ .perm.i, ] <- .result$p.values;
#     verifySizeControl.q.values[ .perm.i, ] <- p.adjust( .result$p.values, method = "fdr" );
#   } # end foreach .perm.i
#   return( list( p.values = verifySizeControl.p.values, q.values = verifySizeControl.q.values ) );
# } # tTestPermPvalWithRefinement.verifySizeControl.454 (..)

#####################chooseSitesToScreen.minimumVariability.454##################################
### TODO: UPDATE/FIX
# chooseSitesToScreen.minimumVariability.454 <- function ( vaccine.aa.freq.list, placebo.aa.freq.list, insert.seq.envir, min.mismatches.to.insert = 4, min.matches.to.insert = min.mismatches.to.insert, min.variant.freq.threshold, normalize.before.filtering = TRUE )
# {
#   insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]];
#   seq.length <- length( insert.char.vector );
#   include.site <- logical( length=seq.length );
#   
#   for( i in 1:seq.length ){
#     vaccine.aa.freq.array <- as.matrix( vaccine.aa.freq.list[[ i ]] );
#     placebo.aa.freq.array <- as.matrix( placebo.aa.freq.list[[ i ]] );
#     all.aa.freq.array <- rbind( vaccine.aa.freq.array, placebo.aa.freq.array );
#     insert.char.pos <- insert.char.vector[ i ];
#     ## TODO: check if this is working properly!!!!!~
#     all.aa.freq.array.normalized <- sapply( 1:nrow( all.aa.freq.array ), 
#                                             function( .row.i ){
#                                               .row.i.aafreq.sum <- sum( all.aa.freq.array[ .row.i, ] );
#                                               if( .row.i.aafreq.sum == 0 ){
#                                                 return( rep( 0, length( all.aa.freq.array[ .row.i, ] ) ) );
#                                               }else{
#                                                 return( all.aa.freq.array[ .row.i, ] / .row.i.aafreq.sum );
#                                               }});
#     
#     # How many subjects have at least one match to the insert?
#     matches.to.insert.pos.by.sequence <- sapply( 1:nrow( all.aa.freq.array.normalized ), function( .sequence ) { all.aa.freq.array.normalized[ .sequence, insert.char.pos ] >= min.variant.freq.threshold } );
#     ### Paul sez: maybe we should instead assert above this that they're not NA?
#     matches.to.insert.pos <- sum( matches.to.insert.pos.by.sequence, na.rm = TRUE );
#     # How many subjects have at least one mismatch to the insert?
#     mismatches.to.insert.pos.by.sequence <- sapply( 1:nrow( all.aa.freq.array.normalized ), function( .sequence ) { any( all.aa.freq.array.normalized[ .sequence, (dimnames( all.aa.freq.array.normalized )[[2]] != insert.char.pos ) ] >= min.variant.freq.threshold ) } );
#     mismatches.to.insert.pos <- sum( mismatches.to.insert.pos.by.sequence, na.rm = TRUE );
#     
#     if( matches.to.insert.pos < min.matches.to.insert ) {
#       ## Not enough insert matches, can't include the position.
#       include.site[ i ] <- FALSE;
#     } else if( mismatches.to.insert.pos < min.mismatches.to.insert ) {
#       ## Not enough insert MISmatches, can't include the position!
#       include.site[ i ] <- FALSE;
#     } else {
#       ## Juuuuust right ;)
#       include.site[ i ] <- TRUE;
#     }
#   }  
#   return( include.site );
# } # chooseSitesToScreen.minimumVariability.454 (..)
