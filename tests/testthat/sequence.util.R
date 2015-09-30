######################## modified from gwj_methods.R ######################
### From "CreateBootstrapDatasets.R":
## Takes a file name of the PAM-style substitution probability matrix, with rows and columns having AAs in alphabetical order (the last row/col is for the gap symbol).  The input file should have no headers and should be delimited by a single space.
## Returns the matrix, with rownames and colnames set appropriately.
readBetty <- function ( .file.name ) {
  .AA.chars <- c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" );
  .gap.chars <- c( "-" );
  .AA.and.gap.chars <- c( .AA.chars, .gap.chars );
  Betty <- read.delim( .file.name, sep=" ", header=F );
  rownames( Betty ) <- .AA.and.gap.chars;
  colnames( Betty ) <- .AA.and.gap.chars;
  ## Sometimes the input file has 0s everywhere for gaps.  If so, make it prob 1 to stay a gap.
  if( Betty[ "-", "-" ] == 0 ) {
    Betty[ "-", "-" ] <- 1;
  }
  return( as.matrix( Betty ) );
} ## readBetty( .. )

#### NOTE HACK! USING STOP CODON SCORE FOR GAP SCORE!
readPAM <- function ( .file.name ) {
  PAM.file <- file( .file.name, "r" );

  have.header <- FALSE;

  PAM.as.list <- list();

  line.text <- readLines( PAM.file, 1 );
  while( length( line.text ) > 0 ) {
    if( line.text == "" ) {
      # blank line.  Do nothing but advance...
      line.text <- readLines( PAM.file, 1 );
      next;
    }
    if( length( grep( "^\\s*#", line.text ) ) > 0 ) {
      # comment line.  Do nothing but advance...
      line.text <- readLines( PAM.file, 1 );
      next;
    }
    if( !have.header ) {
      # Then this should be the header line with the colnames.
      header <- unlist( strsplit( line.text, "\\s+" ) )[ -1 ];
      have.header <- TRUE;
      line.text <- readLines( PAM.file, 1 );
      next;
    }
    line.data <- unlist( strsplit( line.text, "\\s+" ) );
    line.letter <- line.data[ 1 ];
    line.data <- as.numeric( line.data[ -1 ] );
    names( line.data ) <- header;
    PAM.as.list[[ line.letter ]] <- line.data;
    
    # Now prepare for the next round.
    line.text <- readLines( PAM.file, 1 );
  } # End while( length( line.text ) > 0 )
    
  close( PAM.file );

  .AA.chars <- c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" );
  .gap.chars <- c( "-" );
  .AA.and.gap.chars <- c( .AA.chars, .gap.chars );

  ## The PAM header will have a stop codon symbol instead of a gap, but we'll use that for the gap score.
  #### NOTE HACK! USING STOP CODON SCORE FOR GAP SCORE!
  PAM <- matrix( nrow=length( .AA.and.gap.chars ), ncol = length( .AA.and.gap.chars ) );
  rownames( PAM ) <- .AA.and.gap.chars;
  colnames( PAM ) <- .AA.and.gap.chars;
  .stopToGap <- function ( .char ) { if( .char == "*" ) { return(  "-" ); } else { return( .char ); } };
  .gapsToStops <- function ( .vector.of.chars ) { as.vector( sapply( .vector.of.chars, function( .char ) { if( .char == "-" ) { return( "*" ); } else { return( .char ); } } ) ) };
  lapply( names( PAM.as.list ), function( .original.char ) { .char <- .stopToGap( .original.char ); if( .char %in% .AA.and.gap.chars ) { PAM[ .char,  ] <<- PAM.as.list[[ .original.char ]][ .gapsToStops( .AA.and.gap.chars ) ]; } else { NA } } );
  
  return( PAM );
} ## readPAM( .. )


### From Youyi Fong's "youtil.R":
# paste two strings together
# e.g. "a" %+% "b"
"%+%" <- function( a, b )
{
  if( is.numeric( a ) & is.numeric( b ) ) {
    return( sum( a, b ) );
  }
  return( paste( a, b, sep="" ) );
} # End function "%+%"( .. )

### 'aaCountsOnePosition' returns frequencies of 'acceptable.chars' in the character vector 'aa.char.vector'
aaCountsOnePosition <- function ( aa.char.vector, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{
  return( sapply( acceptable.chars, function( .char ){ length( grep( .char, aa.char.vector ) ) } ) );
} ## aaCountsOnePosition(..)

### 'aaCountsPerPosition' returns a list with "aaCountsOnePosition" for each position.
aaCountsPerPosition <- function ( seq.envir, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{ 
  seq.length <- nchar( get(ls(name=seq.envir)[1], env=seq.envir) );
  return( lapply( 1:seq.length, function( .pos ) { return( aaCountsOnePosition( aaPerPosition( seq.envir, .pos ), acceptable.chars=acceptable.chars ) ); } ) );
} # aaCountsPerPosition (..)

### 'aaProportionsPerPosition' is just like aaCountsPerPosition except that it divides each position's table by the total count for that position.
aaProportionsPerPosition <- function ( seq.envir, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{ 
  seq.length <- nchar( get(ls(name=seq.envir)[1], env=seq.envir) );
  aa.proportions <- ( lapply( 1:seq.length, function( .pos ) { .counts <- aaCountsOnePosition( aaPerPosition( seq.envir, .pos ), acceptable.chars=acceptable.chars ); return( .counts / sum( .counts ) ) } ) );
  #num.seqs <- length(list(seq.envir))
  num.seqs <- length(as.list(seq.envir))
  return(list(seq.length = seq.length, aa.proportions = aa.proportions, num.seqs= num.seqs ))
} # aaProportionsPerPosition (..)

### 'aaCountsTwoPositions' returns a two-by-two table of the frequencies of 'acceptable.chars' in the character vector 'aa.char.vector1' and 'aa.char.vector2'
aaCountsTwoPositions <- function ( aa.char.vector1, aa.char.vector2, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ) )
{
  stopifnot( length( aa.char.vector1 ) == length( aa.char.vector2 ) );
  return( sapply( acceptable.chars, function( .char1 ){ sapply( acceptable.chars, function( .char2 ) { length( intersect( grep( .char1, aa.char.vector1 ), grep( .char2, aa.char.vector2 ) ) ) } ) } ) );
} ## aaCountsTwoPositions(..)

# Create the "include.site" vector for use in screenSites(..), based on a minimum-variability criterion (by default, at least 3 sets have to have a seq with an AA the same as the insert char at each included position, and at least 3 sets have to have a seq with an AA that is different from that insert char).  By default, one "set" is just one sequence, but you can provide non-null values of vaccine.sequence.sets.list or placebo.sequence.sets.list if you want to group sequences into sets (eg for multiple sequences from one subject).
# (Optionally, if alwaysScreenInsertGaps) always screen out sites with gap chars ("-") in the insert sequence; this is to acknowledge that the substitution matrices we've been using assign a 0 "probability" score to _any_ substitution from a gap (even to another gap), so it's impossible for GWJ-style methods (t-tests of weights) to find anything there.
chooseSitesToScreen.minimumVariability <- function ( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, min.mismatches.to.insert = 3, min.matches.to.insert = min.mismatches.to.insert, alwaysScreenInsertGaps = FALSE )
{
  insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]];
  seq.length <- length( insert.char.vector );
  include.site <- logical( length=seq.length );

  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.seq.envir ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.seq.envir );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    stopifnot( !is.null( placebo.seq.envir ) );
    placebo.sequence.sets.list = lapply( ls( env=placebo.seq.envir ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.seq.envir );
  }
  
  vaccine.num.sequence.sets <- length( vaccine.sequence.sets.list );
  placebo.num.sequence.sets <- length( placebo.sequence.sets.list );

  for( .pos in 1:seq.length ){
    vaccine.chars.pos <- aaPerPosition( vaccine.seq.envir, .pos );
    placebo.chars.pos <- aaPerPosition( placebo.seq.envir, .pos );
    all.chars.pos <- c( vaccine.chars.pos, placebo.chars.pos );
    insert.char.pos <- insert.char.vector[ .pos ];
    
    if( alwaysScreenInsertGaps && ( insert.char.pos == "-" ) ) {
      include.site[ .pos ] <- FALSE;
    } else {
      # How many subjects have at least one match to the insert?
      matches.to.insert.pos.by.sequence.set <- sapply( c( vaccine.sequence.sets.list, placebo.sequence.sets.list ), function( .sequences.in.set ) { as.numeric( any( all.chars.pos[ .sequences.in.set ] == insert.char.pos, na.rm = TRUE ) ) } );
      matches.to.insert.pos <- sum( matches.to.insert.pos.by.sequence.set, na.rm=T );
      # How many subjects have at least one mismatch to the insert?
      mismatches.to.insert.pos.by.sequence.set <- sapply( c( vaccine.sequence.sets.list, placebo.sequence.sets.list ), function( .sequences.in.set ) { as.numeric( any( all.chars.pos[ .sequences.in.set ] != insert.char.pos, na.rm = TRUE ) ) } );
      mismatches.to.insert.pos <- sum( mismatches.to.insert.pos.by.sequence.set, na.rm=T );
      
      if( matches.to.insert.pos < min.matches.to.insert ) {
        include.site[ .pos ] <- FALSE;
      } else if( mismatches.to.insert.pos < min.mismatches.to.insert ) {
        include.site[ .pos ] <- FALSE;
      } else {
        include.site[ .pos ] <- TRUE;
      }
    } # End if alwaysScreenInsertGaps && ( insert.char.pos == "-" ) .. else ..
  }
  
  return( include.site );
} # chooseSitesToScreen.minimumVariability (..)

### 'screenOutSites' returns an environment with sequences containing sites specified in 'include.site'
# See also updateSiteSetsList(..)
screenOutSites <- function(seq.envir, include.site)
{
  stopifnot( is.environment( seq.envir ) );
	seq.length <- nchar(get(ls(name=seq.envir)[1], envir=seq.envir))
  if (seq.length != length(include.site)){ stop(paste( "Sequence length and 'include.site' vector length differ: seq.length=", seq.length, ", length(include.site)=", length(include.site), ".", sep="" )) }
	for (seqid in ls(name=seq.envir)){
		seq <- get(seqid, env=seq.envir)
		seq <- paste(strsplit(seq,"")[[1]][include.site], collapse="")
		assign(seqid, seq, env=seq.envir)
	}
	return (seq.envir )
} ## screenOutSites(..)

# Since screenOutSites(..) actually removes sites from sequences, we might lose track of our numbering scheme.  Eg what was site 50 before screening out sites 1-49 will become site 1 afterwards.  This function takes a list of site-vectors, with each site falling in the range 1..length( include.site ).  It returns the same list, but with the site values updated to the new numbering scheme after removing all sites for which include.site is FALSE.  Note that sites that are not included will become NAs in the result list, so you might have site sets that are entirely NAs.
updateSiteSetsList <- function ( include.site, site.sets.list = NULL )
{
  if( is.null( site.sets.list ) ) {
    site.sets.list <- lapply( which( include.site ), function( .site ) { .site } );
  }
  old.site.to.new.site <- cumsum( as.numeric( include.site ) );
  old.site.to.new.site[ !include.site ] <- NA;
  return( lapply( site.sets.list, function( .site.set ) { old.site.to.new.site[ .site.set ] } ) );
} # updateSiteSetsList (..)

##########
### From "CreateBootstrapDatasets.R":
readFastaFile <- function ( .file.name, .acceptable.characters = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), .replace.unacceptable.characters = '-', .parent = parent.frame(), .seqid.sep = " ", .seqid.length.limit = 1000, be.verbose = TRUE )
{
  if( be.verbose ) {
    cat( "Reading fasta file: " %+% .file.name, fill=TRUE );
  }
  
  return.hash <- new.env( hash=TRUE, parent=.parent );
  .fastaFile <- file( .file.name, "r" );

  .acceptable.characters.asonestring <- toupper( paste( .acceptable.characters, sep="", collapse="" ) );
  
  .line.text <- readLines( .fastaFile, 1 );
  .seqid <- NULL;
  while( length( .line.text ) > 0 ) { # Outer loop looks for new sequence blocks (starting w a seqid)
    if( !grep( "\\S", .line.text ) ) { # Just whitespace
      ## Skip whitespace lines.
      .line.text <- readLines( .fastaFile, 1 );
      next;
    }
    stopifnot( substr( .line.text, 1, 1 ) == ">" ); ## We expect sequence blocks to begin with a descripotion line
    .seqid <- gsub( ">\\s*(\\S.*)", "\\1", .line.text ); ## Ignore ws between ">" and seqid.
    .seqid <- substr( strsplit( .seqid, .seqid.sep )[[1]][1], 1, .seqid.length.limit ); # Take first .seqid.length.limit chars up to first .seqid.sep char
    if( exists( .seqid, inherits=TRUE, envir=return.hash  ) ) { # If the .seqid already has an associated sequence
      warning( "Duplicate seqids (" %+% .seqid %+% ").  Only the first instance will be used!" );
      .line.text <- readLines( .fastaFile, 1 );
      next;
    }
    .line.text <- readLines( .fastaFile, 1 );
    .temp.seq <- "";
    while( length( .line.text ) > 0 ) { # Inner loop reads sequence associated with this block
      if( substr( .line.text, 0, 1 ) == ">" ) { # Found the next block / seqid.
        break; # Break from inner loop.
      }
      .temp.seq <- .temp.seq %+% .line.text;
      .line.text <- readLines( .fastaFile, 1 );
    } # End inner loop (reading in the lines of the current sequence)
    if( is.na( .replace.unacceptable.characters ) ) {
      stopifnot( grep( "[^" %+% .acceptable.characters.asonestring %+% "]", .temp.seq, invert=T ) );
    } else {
      .temp.seq <- gsub( "[^" %+% .acceptable.characters.asonestring %+% "]", .replace.unacceptable.characters, .temp.seq );
    }
    assign( .seqid, .temp.seq, envir=return.hash );
  } # End while( length( .line.text ) > 0 )

  close( .fastaFile );
  return( return.hash );
} # readFastaFile (..)

## Modified from Youyi Fong's youtil.R
## Takes a list that maps seqid (strings) to aligned sequence (strings) and a filename, and optionally a prefix to be appended to the output sequences, and an optional first sequence number (if non-NA, a number will be prepended to the prefix and incremented for each sequence, to help ensure unique seqids).
## Returns nothing.  Writes Fasta-formated sequences as a side effect.
## 
writeFastaFile <- function( .seq.list, .file.name, .seqid.prefix = "", .first.seqid.number = NA )
{
  .outFile <- file( .file.name, open="w" );
  .seq.names <- names( .seq.list );
  .seqid.number <- .first.seqid.number;
  for( .seq.i in 1:length( .seq.list ) ) {
    write( file=.outFile, ">" %+% ifelse( is.na( .seqid.number ), '', as.character( .seqid.number ) %+% "|" ) %+% .seqid.prefix %+% .seq.names[ .seq.i ], append=T );
    if( !is.na( .seqid.number ) ) {
      .seqid.number <- .seqid.number + 1;
    }
    write( file=.outFile, .seq.list[[.seq.i]], append=T );
  }
  close( .outFile );
} # writeFastaFile(..)

### From "CreateBootstrapDatasets.R":
countFastaSeqs <- function ( seq.environment )
{
  return( length( as.list( seq.environment ) ) );
} ## countFastaSeqs (..)

### 'aaPerPosition' returns a character vector of aas at the position 'position' for sequences in 'seq.environment'
## NOTE it'll return any char, including "unacceptable" chars.
aaPerPosition <- function ( seq.environment, position )
{
  #OLD:return( do.call( "c", lapply(as.list(seq.environment), function(seq){ strsplit(seq,"")[[1]][position] }) ) );
  # FASTER:
  #return( unlist( lapply(as.list(seq.environment), function(seq){ strsplit(seq,"")[[1]][position] }) ) );
  # FASTEST:
  return( unlist( sapply( as.list( seq.environment ), function( .seq ) { substr( .seq, position, position ) } ) ) );
  # VARIOUS OTHER WAYS.  Bizarrely, it appears that unlist(lapply(as.list(..),..)) is noticably faster than sapply( as.list(..), .. ).  Unsurprisingly, substr(..) is faster than strsplit(..).
  # return( unlistsapply( mget( ls( env=seq.environment ), env=seq.environment ), function( .seq ) { substr( .seq, position, position ) } ) );
  #return( sapply(as.list(seq.environment), function(seq){ strsplit(seq,"")[[1]][position] }) );
  # return( sapply(mget(ls( env=seq.environment ), env=pla.sample.env), function(seq){ strsplit(seq,"")[[1]][position] }) );
  # return( sapply( mget( ls( env=seq.environment ), env=seq.environment ), function( .seq ) { substr( .seq, position, position ) } ) );
} ## aaPerPosition (..)

### charsToCategories returns a numeric vector of the indices of the given character vector into the given acceptable.chars vector.  Unacceptable chars will be replaced by the given unacceptable.char.value (defaults to length( acceptable.chars ) + 1; meaning that the actual number of categories may be greater than the number of acceptable chars!
charsToCategories <- function ( char.vec, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), unacceptable.char.value = ( length( acceptable.chars ) + 1 ) )
{
  category.vec <- as.numeric( addNA( factor( char.vec, levels=acceptable.chars ) ) );
  if( is.na( unacceptable.char.value ) ) {
    unacceptable.char.value <- length( acceptable.chars ) + 1;
  }
  if( ( unacceptable.char.value != ( length( acceptable.chars ) + 1 ) ) && ( sum( category.vec == ( length( acceptable.chars ) + 1 ) ) > 0 ) ) {
    category.vec[ category.vec == ( length( acceptable.chars ) + 1 ) ] <- unacceptable.char.value;
  }
  return( category.vec );
} ## charsToCategories (..)


filterCharVector <- function ( aa.char.vector, acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" ), .grep.pattern = paste( "[", paste( acceptable.chars, collapse="|" ), "]", sep="" ) )
{
  return( aa.char.vector[ grep( .grep.pattern, aa.char.vector ) ] );
} # filterCharVector (..)

######################## END modified from gwj_methods.R ##################


# This function takes a vector of values and returns a table with counts per categories, but only counts values in the range 1:num.categories.  Any value that isn't an integer or isn't in that range will be ignored.  Note that it's ok to pass NULL for categorical.data.vector (it'll just be a table of 0s).
makeCountsTable <- function ( categorical.data.vector, num.categories )
{
  counts.table <- rep( 0, num.categories );
  names( counts.table ) <- 1:num.categories;
  .tmp <- table( categorical.data.vector );
  .tmp <- .tmp[ intersect( names( counts.table ), names( .tmp ) ) ];
  counts.table[ names( .tmp ) ] <- .tmp;

  return( counts.table );
} # makeCountsTable (..)

##  TODO: Come up with a better name for this one.  Make a pretty counts table for diagnosing sites-found-to-have-sieve-effects.
make.counts.table <- function( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, pos, remove.empty.columns = TRUE, sort.insert.on.left = TRUE, sort.by.pla = TRUE ) {
  .table <- rbind( t( as.matrix( aaCountsOnePosition( aaPerPosition( vaccine.seq.envir, pos ) ) ) ), t( as.matrix( aaCountsOnePosition( aaPerPosition( placebo.seq.envir, pos ) ) ) ) );
  rownames( .table ) <- c( "V", "P" );
  .colnames <- tolower( colnames( .table ) );
  names( .colnames ) <- .colnames;
  insert.char <- substr( get( ls( envir=insert.seq.envir )[ 1 ], envir=insert.seq.envir ), pos, pos );
  .colnames[ tolower( insert.char ) ] <- toupper( insert.char );
  colnames( .table ) <- .colnames;
  if( remove.empty.columns ) {
    .table <- .table[ , apply( .table, 2, sum ) != 0, drop = F ];
    .colnames <- colnames( .table );
  } # End if remove.empty.columns
  if( sort.by.pla && ( length( .colnames ) > 1 ) ) {
    #.table <- .table[ , sort( .table[ "P", ], index.return = T, decreasing = T )$ix ];
    .table <- .table[ , order( .table[ "P", ], .table[ "V", ], decreasing = T ) ];
    .colnames <- colnames( .table );
  }
  if( sort.insert.on.left ) {
    .insert.column <- which( .colnames == toupper( insert.char ) );
    if( length( .insert.column ) == 0 ) {
      # Insert not observed at all
      .table <- cbind( c( 0, 0 ), .table );
      .colnames <- c( insert.char, .colnames );
      colnames( .table ) <- .colnames;
    } else {
      if( .insert.column != 1 ) {
        .table <- cbind( .table[ , .insert.column ], .table[ , -.insert.column ] );
        .colnames <- c( .colnames[ .insert.column ], .colnames[ -.insert.column ] );
        colnames( .table ) <- .colnames;
      }
    }
  } # End if sort.insert.on.left
  return( .table );
} # make.counts.table ( ... )

# pos1 in rows, pos2 in columns
make.counts.covariation.table <- function( breakthrough.seq.envir, insert.seq.envir, pos1, pos2, remove.empty.columns.and.rows = TRUE, sort.insert.on.left.and.top = TRUE ) {
  .covariation.table <- aaCountsTwoPositions( aaPerPosition( breakthrough.seq.envir, pos2 ), aaPerPosition( breakthrough.seq.envir, pos1 ) );
  .colnames <- tolower( colnames( .covariation.table ) );
  names( .colnames ) <- .colnames;
  .rownames <- tolower( rownames( .covariation.table ) );
  names( .rownames ) <- .rownames;
  insert.char1 <- substr( get( ls( envir=insert.seq.envir )[ 1 ], envir=insert.seq.envir ), pos1, pos1 );
  insert.char2 <- substr( get( ls( envir=insert.seq.envir )[ 1 ], envir=insert.seq.envir ), pos2, pos2 );
  .colnames[ tolower( insert.char2 ) ] <- toupper( insert.char2 );
  colnames( .covariation.table ) <- .colnames;
  .rownames[ tolower( insert.char1 ) ] <- toupper( insert.char1 );
  rownames( .covariation.table ) <- .rownames;
  if( remove.empty.columns.and.rows ) {
    .covariation.table <- .covariation.table[ apply( .covariation.table, 1, sum ) != 0, apply( .covariation.table, 2, sum ) != 0, drop = F ];
    .rownames <- rownames( .covariation.table );
    .colnames <- colnames( .covariation.table );
  } # End if remove.empty.columns
  if( sort.insert.on.left.and.top ) {
    .insert.column <- which( .colnames == toupper( insert.char2 ) );
    .insert.row <- which( .rownames == toupper( insert.char1 ) );
    if( length( .insert.column ) == 0 ) {
      # Insert not observed at all in grp 2
      .covariation.table <- cbind( rep( 0, nrow( .covariation.table ) ), .covariation.table );
      .colnames <- c( insert.char2, .colnames );
      colnames( .covariation.table ) <- .colnames;
    } else {
      # Rearrange so the insert column is on the left.
      if( .insert.column != 1 ) {
        .covariation.table <- cbind( .covariation.table[ , .insert.column, drop = FALSE ], .covariation.table[ , -.insert.column, drop = FALSE ] );
        .colnames <- c( .colnames[ .insert.column ], .colnames[ -.insert.column ] );
        colnames( .covariation.table ) <- .colnames;
      }
    } # End if length( .insert.column ) == 0 .. else ..
    if( length( .insert.row ) == 0 ) {
      # Insert not observed at all in grp 1
      .covariation.table <- rbind( rep( 0, ncol( .covariation.table ) ), .covariation.table );
      .rownames <- c( insert.char1, .rownames );
      rownames( .covariation.table ) <- .rownames;
    } else {
      # Rearrange so the insert row is on the top.
      if( .insert.row != 1 ) {
        .covariation.table <- rbind( .covariation.table[ .insert.row, , drop = FALSE ], .covariation.table[ -.insert.row, , drop = FALSE ] );
        .rownames <- c( .rownames[ .insert.row ], .rownames[ -.insert.row ] );
        rownames( .covariation.table ) <- .rownames;
      }
    } # End if length( .insert.row ) == 0 .. else ..
  } # End if sort.insert.on.left.and.top
  if( ( length( .colnames ) > 1 ) & ( length( .rownames ) > 1 ) ) {
    if( sort.insert.on.left.and.top ) {
      .covariation.table <- .covariation.table[ c( 1, 1 + order( .covariation.table[ -1, 1 ], .covariation.table[ -1, 2 ], decreasing = T ) ), c( 1, 1 + order( .covariation.table[ 1, -1 ], .covariation.table[ 2, -1 ], decreasing = T ) ) ];
    } else {
      .covariation.table <- .covariation.table[ order( .covariation.table[ , 1 ], .covariation.table[ , 2 ], decreasing = T ), order( .covariation.table[ 1, ], .covariation.table[ 2, ], decreasing = T ) ];
    }
  }
  return( .covariation.table );
} # make.counts.covariation.table ( ... )

make.counts.differential.covariation.table <- function( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, pos1, pos2, remove.empty.columns.and.rows = TRUE, sort.insert.on.left.and.top = TRUE ) {
  return( list( vaccine = make.counts.covariation.table( vaccine.seq.envir, insert.seq.envir, pos1, pos2, remove.empty.columns.and.rows=remove.empty.columns.and.rows, sort.insert.on.left.and.top=sort.insert.on.left.and.top ), placebo = make.counts.covariation.table( placebo.seq.envir, insert.seq.envir, pos1, pos2, remove.empty.columns.and.rows=remove.empty.columns.and.rows, sort.insert.on.left.and.top=sort.insert.on.left.and.top ) ) );
} # make.counts.differential.covariation.table ( ... )

fractions.table.from.counts.table <- function ( counts.table )
{
  counts.table / rowSums( counts.table )
} # fractions.table.from.counts.table (..)

plot.fractions.table <- function ( fractions.table, colors.for.AAs.fn, ... )
{
  barplot( t( fractions.table ), col = colors.for.AAs.fn( toupper( colnames( fractions.table ) ) ), ... );
} # plot.fractions.table (..)

pretty.plot.fractions.table <- function ( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, pos, include.legend = FALSE, exclude.unrepresented.properties.from.legend = TRUE, sort.by.color = FALSE, label = NA, sublabel = NA, color.scheme = "fusheng" )
{
  ## Set up colors
    
  #MAEditor.aa.colors <- list( A = "light green", G = "light green", C = "green", D = "dark green", E = "dark green", N = "dark green", Q = "dark green", I = "blue", L = "blue", M = "blue", V = "blue", F = "violet", W = "violet", Y = "violet", H = "dark blue", K = "orange", R = "orange", P = "pink", S = "red", T = "red" );
  #cinema.aa.colors.by.property <- list( polar.positive.HKR. = "cyan", polar.negative.DE. = "orange", polar.neutral.STNQ. = "green", non.polar.aliphatic.AVLIM. = "white", non.polar.aromatic.FWY. = "red", non.polar.aromatic.PG. = "brown", non.polar.aromatic.C. = "pink", special = "grey" );
  #cinema.aa.property.map <- list( H = "polar.positive.HKR.", K = "polar.positive.HKR.", R = "polar.positive.HKR.", D = "polar.negative.DE.", E = "polar.negative.DE.", S = "polar.neutral.STNQ.", T = "polar.neutral.STNQ.", N = "polar.neutral.STNQ.", Q = "polar.neutral.STNQ.", A = "non.polar.aliphatic.AVLIM.", V = "non.polar.aliphatic.AVLIM.", L = "non.polar.aliphatic.AVLIM.", I = "non.polar.aliphatic.AVLIM.", M = "non.polar.aliphatic.AVLIM.", F = "non.polar.aromatic.FWY.", W = "non.polar.aromatic.FWY.", Y = "non.polar.aromatic.FWY.", P = "non.polar.aromatic.PG.", G = "non.polar.aromatic.PG.", C = "non.polar.aromatic.C.", X = "special" )
  #cinema.aa.property.map[[ "-" ]] <- "special";
  #cinema.aa.colors <- lapply( cinema.aa.property.map, function( .property ) { cinema.aa.colors.by.property[[ .property ]] } )

  if( color.scheme == "fusheng" ) {
    fusheng.aa.colors.by.property <- list( sulfhydryl.C. = "orange", small.hydrophilic.STPAG. = "white", acid.NDEQ. = "red", basic.HRK. = "cyan", small.hydrophobic.MILV. = "yellow", aromatic.FYW. = "pink", special = "grey" )
  #   sulfhydryl (C),
  # small hydrophilic (STPAG),
  # acid (NKEQ), --> K -> D --> NDEQ
  # basic (HRK),
  # small hydrophobic (MILV)
  # aromatic (FYW).
    fusheng.aa.property.map <- list( H = "basic.HRK.", R = "basic.HRK.", K = "basic.HRK.", N = "acid.NDEQ.", D = "acid.NDEQ.", E = "acid.NDEQ.", Q = "acid.NDEQ.", S = "small.hydrophilic.STPAG.", T = "small.hydrophilic.STPAG.", P = "small.hydrophilic.STPAG.", A = "small.hydrophilic.STPAG.", G = "small.hydrophilic.STPAG.", M = "small.hydrophobic.MILV.", I = "small.hydrophobic.MILV.", L = "small.hydrophobic.MILV.", V = "small.hydrophobic.MILV.", F = "aromatic.FYW.", Y = "aromatic.FYW.", W = "aromatic.FYW.", C = "sulfhydryl.C.", X = "special" );
    fusheng.aa.property.map[[ "-" ]] <- "special";
    fusheng.aa.colors <- lapply( fusheng.aa.property.map, function( .property ) { fusheng.aa.colors.by.property[[ .property ]] } )
  
    aa.colors.by.property <- fusheng.aa.colors.by.property;
    aa.property.map <- fusheng.aa.property.map;
    aa.colors <- fusheng.aa.colors;
  } else if( color.scheme == "rainbow" ) {
    .AA.chars <- c( "H", "R", "K", "N", "D", "E", "Q", "S", "T", "P", "A", "G", "M", "I", "L", "V", "F", "Y", "W", "C" );
    .gap.chars <- c( "-" );
    .AA.and.gap.chars <- c( .AA.chars, .gap.chars );

    rainbow.aa.colors <- c( rainbow( 20 ), "grey" );
    names( rainbow.aa.colors ) <- .AA.and.gap.chars;
    rainbow.aa.colors.by.property <- rainbow.aa.colors;
    rainbow.aa.property.map <- .AA.and.gap.chars;
    names( rainbow.aa.property.map ) <- .AA.and.gap.chars;
    
    aa.colors.by.property <- rainbow.aa.colors.by.property;
    aa.property.map <- rainbow.aa.property.map;
    aa.colors <- rainbow.aa.colors;
  } else {
     # If not "fusheng" or "rainbow" use the color.scheme as the color for all colors.
     # Or if unspecified, use "white"
    if( is.na( color.scheme ) || ( color.scheme == "" ) ) {
      color.scheme <- "white";
    }
    .AA.chars <- c( "H", "R", "K", "N", "D", "E", "Q", "S", "T", "P", "A", "G", "M", "I", "L", "V", "F", "Y", "W", "C" );
    .gap.chars <- c( "-" );
    .AA.and.gap.chars <- c( .AA.chars, .gap.chars );

    default.aa.colors <- rep( color.scheme, 21 );
    names( default.aa.colors ) <- .AA.and.gap.chars;
    default.aa.colors.by.property <- default.aa.colors;
    default.aa.property.map <- .AA.and.gap.chars;
    names( default.aa.property.map ) <- .AA.and.gap.chars;
    
    aa.colors.by.property <- default.aa.colors.by.property;
    aa.property.map <- default.aa.property.map;
    aa.colors <- default.aa.colors;
  }
  
  .colors.for.AAs.fn <- function( .vector.of.AAs ) { unlist( lapply( .vector.of.AAs, function( .aa ) { aa.colors[[ toupper( .aa ) ]] } ) ) }
  
  .counts.table <- make.counts.table( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, pos );

  # Sort by color?  Insert color will be on the left; sort colors by total count (in pla) and within a color, keep sort order as before [which is sorted by count (in pla)].
  if( sort.by.color ) {
    # Gather AAs with the same colors
    # On the left, use the insert AA; then all AAs of the same color.
    .AAs <- colnames( .counts.table );
    .colors.for.AAs <- .colors.for.AAs.fn( .AAs );
    names( .colors.for.AAs ) <- .AAs;
    .insert.AA <- .AAs[ 1 ];
    .insert.color <- .colors.for.AAs[ 1 ];
    .AA.indices.with.same.color.as.insert <- which( .colors.for.AAs == .insert.color );
    .new.AAs <- .AAs[ .AA.indices.with.same.color.as.insert ];
    # We've put those on the left.  Now group the rest by color.
    .rest <- .AAs[ -.AA.indices.with.same.color.as.insert ];
    if( length( .rest ) > 0 ) {
      .rest.colors.sorted <- sort( .colors.for.AAs[ .rest ] );
      # Sort those groups-of-AAs-with-same-colors by total count in pla group.
      .total.count.for.color.pla <- sort( sapply( unique( .rest.colors.sorted ), function( .color ) { sum( .counts.table[ "P", .AAs[ .colors.for.AAs == .color ] ] ) } ), decreasing = TRUE );
      for( .unique.color.i in 1:length( .total.count.for.color.pla ) ) {
        .new.AAs <- c( .new.AAs, names( .rest.colors.sorted )[ .rest.colors.sorted == names( .total.count.for.color.pla )[ .unique.color.i ] ] );
      }
    }
    .counts.table <- .counts.table[ , .new.AAs, drop = FALSE ];
  } # End if sort.by.color

  .fractions.table <- fractions.table.from.counts.table( .counts.table )
  .v.c <- .counts.table[ "V", ]
  .v.c <- .v.c[ .v.c != 0 ]
  .p.c <- .counts.table[ "P", ]
  .p.c <- .p.c[ .p.c != 0 ]
  
  .x.locations.of.bars <- plot.fractions.table( .fractions.table, .colors.for.AAs.fn, ylab = "fraction of subjects", width = rep( .1, 2 ), xlim = c( 0, 1 ), ylim = c( -0.1, 1.1 ), axisnames = F )
  text( .x.locations.of.bars, -0.05, c( "V", "P" ), srt = 90, font = 2 )

  if( !is.na( label ) ) {
    text( mean( .x.locations.of.bars ), 1.08, label, font = 2 )
  }
  if( !is.na( sublabel ) ) {
    text( mean( .x.locations.of.bars ), 1.03, sublabel, font = 1 )
  }

  if( length( .v.c ) > 1 ) {
    text( 0, c( 0, ( cumsum( .v.c ) / sum( .v.c ) )[ -length( .v.c ) ] ), toupper( names( .v.c ) ), adj = c( .05, .75 ), srt = 90, col = "black", cex = .6, font = 2 )
    .v.c.noones <- as.character( .v.c )
    .v.c.noones[ .v.c == 1 ] <- ""
    text( .1, c( 0, ( cumsum( .v.c ) / sum( .v.c ) )[ -length( .v.c ) ] ), .v.c.noones, adj = c( -.15, .5 ), srt = 90, col = "black", cex = 1.0, font = 2 )
  }
  if( length( .p.c ) > 1 ) {
    text( .x.locations.of.bars[ 2 ] + .06, c( 0, ( cumsum( .p.c ) / sum( .p.c ) )[ -length( .p.c ) ] ), toupper( names( .p.c ) ), adj = c( .05, 1.0 ), srt = 90, col = "black", cex = .6, font = 2 )
    .p.c.noones <- as.character( .p.c )
    .p.c.noones[ .p.c == 1 ] <- ""
    text( .x.locations.of.bars[ 2 ] - .04, c( 0, ( cumsum( .p.c ) / sum( .p.c ) )[ -length( .p.c ) ] ), .p.c.noones, adj = c( -.15, 1 ), srt = 90, col = "black", cex = 1.0, font = 2 )
  }

  ## Make the legend
  if( include.legend ) {
    if( exclude.unrepresented.properties.from.legend ) {
      .AAs <- colnames( .counts.table );
      aa.properties <- unique( aa.property.map[ toupper( .AAs ) ] );
    } else {
      aa.properties <- names( aa.colors.by.property );
    }
    aa.properties.txt <- gsub( "non\\.", "non-", aa.properties );
    aa.properties.txt <- gsub( "\\.([^\\.]+)\\.$", " (\\1)", aa.properties.txt );
    aa.properties.txt <- gsub( "\\.", " ", aa.properties.txt );
    aa.properties.txt <- gsub( "special", "gap or missing data", aa.properties.txt );
    .height <- .2;
    for( .aa.properties.i in 1:length( aa.properties ) ) {
      .last.result <- legend( .5, .height, aa.properties.txt[ .aa.properties.i ], text.col = "black", bg = aa.colors.by.property[[ aa.properties[[ .aa.properties.i ]] ]], text.width = .38 );
      .height <- .height + .last.result$rect$h;
    }
  } # End if include.legend 
} # pretty.plot.fractions.table (..)

hxb2Pos.to.alignmentPos <- function ( hxb2Pos, hxb2.map )
{
  hxb2.map[ hxb2.map[ , "hxb2Pos" ] %in% pos, "posNum" ]
}

alignmentPos.to.hxb2Pos <- function ( alignmentPos, hxb2.map )
{
  as.character( hxb2.map[ alignmentPos, "hxb2Pos" ] )
}

############ Via Craig Magaret ########
returnLBDonorVectorDict <- function ()
{
   # we need a dictionary to define the LB vector values; let's store it
   # here centrally -- to be called by subroutine -- so that we don't have
   # multiple copies floating around
   .AA.chars <- c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" );
  .gap.chars <- c( "-" );
  .AA.and.gap.chars <- c( .AA.chars, .gap.chars );

  LBDonorVectorDict <- matrix( nrow = 21, ncol = 12 );
  rownames( LBDonorVectorDict ) <- .AA.and.gap.chars;
  colnames( LBDonorVectorDict ) <- c( 
    "Hydrophobic",
    "Polar",
    "Small",
    "Proline",
    "Tiny",
    "Aliphatic",
    "Aromatic",
    "Positive",
    "Negative",
    "Charged",
    "Hydrogen Donors",
    "Hydrogen Acceptors"
  );
  LBDonorVectorDict[ "-", ] <- c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 );
  LBDonorVectorDict[ "A", ] <- c( 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 );
  #LBDonorVectorDict[ "B", ] <- c( 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "C", ] <- c( 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "D", ] <- c( 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1 );
  LBDonorVectorDict[ "E", ] <- c( 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1 );
  LBDonorVectorDict[ "F", ] <- c( 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "G", ] <- c( 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "H", ] <- c( 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1 );
  LBDonorVectorDict[ "I", ] <- c( 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "K", ] <- c( 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0 );
  LBDonorVectorDict[ "L", ] <- c( 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "M", ] <- c( 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "N", ] <- c( 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "P", ] <- c( 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "Q", ] <- c( 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "R", ] <- c( 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0 );
  LBDonorVectorDict[ "S", ] <- c( 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "T", ] <- c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "V", ] <- c( 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 );
  LBDonorVectorDict[ "W", ] <- c( 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1 );
  LBDonorVectorDict[ "Y", ] <- c( 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1 );
  #LBDonorVectorDict[ "Z", ] <- c( 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 );

   return( LBDonorVectorDict );
} # returnLBDonorVectorDict (..)

# Return a 0-1 weight matrix with identical rows (that is, not insert-dependent) with 0 meaning that the column-AA has the property, and 1 meaning that it does not.  Note that the gap column will have a "1" but the gap row will be just like the other rows.
makeSimilarityMatrixFromLBPhysicochemicalProperty <-
  function ( property, LBDonorVectorDict = returnLBDonorVectorDict() )
{
  .LBProperties <- colnames( LBDonorVectorDict );
  .AA.and.gap.chars <- rownames( LBDonorVectorDict );
  stopifnot( property %in% .LBProperties );
  .matrix <- rbind( matrix( rep( LBDonorVectorDict[ , property ], times = 21 ), nrow = 21, byrow = TRUE ) );
  rownames( .matrix ) <- .AA.and.gap.chars;
  colnames( .matrix ) <- .AA.and.gap.chars;
  # Make sure all gap weights are 0.
  .matrix[ , "-" ] <- 0;
  return( .matrix );
} # makeSimilarityMatrixFromLBPhysicochemicalProperty (..)

# Note: this returns a symmetric weight matrix with non-negative integer values.
makeWeightMatrixCountingSharedLBPhysicochemicalProperties <- function ( LBDonorVectorDict = returnLBDonorVectorDict() )
{
  WeightMatrix <- matrix( nrow = nrow( LBDonorVectorDict ), ncol = nrow( LBDonorVectorDict ) );
  rownames( WeightMatrix ) <- rownames( LBDonorVectorDict );
  colnames( WeightMatrix ) <- rownames( LBDonorVectorDict );
  for( from.i in 1:nrow( WeightMatrix ) ) {
    for( to.i in from.i:nrow( WeightMatrix ) ) {
      properties.in.common <- sum( LBDonorVectorDict[ from.i, ] & LBDonorVectorDict[ to.i, ] );
      WeightMatrix[ from.i, to.i ] <- properties.in.common;
      WeightMatrix[ to.i, from.i ] <- properties.in.common;
    } # End foreach to.i
  } # End foreach from.i
  return( WeightMatrix );
} # makeWeightMatrixCountingSharedLBPhysicochemicalProperties (..)

makeProbabilityMatrixProportionalToCountsMatrix <- function( CountsMatrix, minimum.probability = 1E-5 )
{
  # First convert any 0s to NA, then normalize rest to sum to 1-num.NAs
  return( t( apply( CountsMatrix, 1, function( .row ) { .row[ .row == 0 ] <- NA; .num.NAs <- sum( is.na( .row ) ); .row <- ( .row / sum( .row, na.rm = TRUE ) ) * ( 1 - ( .num.NAs * minimum.probability ) ); .row[ is.na( .row ) ] <- minimum.probability; return( .row ); } ) ) );
} # makeProbabilityMatrixProportionalToCountsMatrix (..)


## Returns vector of the #sequences at each site that have no matches to the vaccine insert.
sitesWithNoMatchToInsert <- function ( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL )
{
  insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]];
  seq.length <- length( insert.char.vector );
  num.sequences.with.no.match <- rep( NA, length=seq.length );

    if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.seq.envir ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.seq.envir );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    stopifnot( !is.null( placebo.seq.envir ) );
    placebo.sequence.sets.list = lapply( ls( env=placebo.seq.envir ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.seq.envir );
  }
  
  vaccine.num.sequence.sets <- length( vaccine.sequence.sets.list );
  placebo.num.sequence.sets <- length( placebo.sequence.sets.list );

  for( i in 1:seq.length ){
    vaccine.chars.pos <- aaPerPosition( vaccine.seq.envir, i );
    placebo.chars.pos <- aaPerPosition( placebo.seq.envir, i );
    all.chars.pos <- c( vaccine.chars.pos, placebo.chars.pos );
    insert.char.pos <- insert.char.vector[ i ];

    # Which sites have no sequences that match the vaccine insert?
    num.sequences.with.no.match[ i ] <- sum( !(all.chars.pos == insert.char.pos), na.rm = TRUE );                                 
  }
  return( num.sequences.with.no.match );
} # sitesWithNoMatchToInsert (..)
  
# TODO: rename to be more descriptive
# Returns the #(included sites) after calling chooseSitesToScreen.minimumVariability (..) that have a minimum number of noninsert AAs (and insert AAs) observed at each site.
numIncludedSitesAfterMinimumVariabilityFilter <- function(  min.num.noninsert.aas = 3, vaccine.fasta.file.name, placebo.fasta.file.name, insert.fasta.file.name, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL, be.verbose = TRUE )
  {      
    if( is.environment( insert.fasta.file.name ) ) {
      insert.env <- insert.fasta.file.name;
    } else {
      insert.env <- readFastaFile( insert.fasta.file.name, be.verbose=be.verbose );
    }
    if( is.environment( vaccine.fasta.file.name ) ) {
      vaccine.env <- vaccine.fasta.file.name;
    } else {
      vaccine.env <- readFastaFile( vaccine.fasta.file.name, be.verbose=be.verbose );
    }
    if( is.environment( placebo.fasta.file.name ) ) {
      placebo.env <- placebo.fasta.file.name;
    } else {
      placebo.env <- readFastaFile( placebo.fasta.file.name, be.verbose=be.verbose );
    }
    vaccine.sequence.sets.list <- NULL;
    placebo.sequence.sets.list <- NULL;
    
    min.var.filter.include.site <- chooseSitesToScreen.minimumVariability( vaccine.env, placebo.env, insert.env, vaccine.sequence.sets.list=vaccine.sequence.sets.list, placebo.sequence.sets.list=placebo.sequence.sets.list, min.mismatches.to.insert = min.num.noninsert.aas );
    num.included.sites <- sum( min.var.filter.include.site, na.rm = TRUE );
    return( num.included.sites );
  }

## returns a table of the number of sequences that match/mismatch the insert AA, separated by vaccine/placebo group.
makeMatchMismatchInsertAAtable <- function ( vaccine.seq.envir, placebo.seq.envir, insert.seq.envir, vaccine.sequence.sets.list = NULL, placebo.sequence.sets.list = NULL )
{
  insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]];
  seq.length <- length( insert.char.vector );
  match.mismatch.table <- array(NA, dim=c( seq.length, 6 ));

  if( is.null( vaccine.sequence.sets.list ) ) {
    vaccine.sequence.sets.list = lapply( ls( env=vaccine.seq.envir ), function( .seq.name ) { .seq.name } );
    names( vaccine.sequence.sets.list ) <- ls( env=vaccine.seq.envir );
  }
  if( is.null( placebo.sequence.sets.list ) ) {
    stopifnot( !is.null( placebo.seq.envir ) );
    placebo.sequence.sets.list = lapply( ls( env=placebo.seq.envir ), function( .seq.name ) { .seq.name } );
    names( placebo.sequence.sets.list ) <- ls( env=placebo.seq.envir );
  }
  
  vaccine.num.sequence.sets <- length( vaccine.sequence.sets.list );
  placebo.num.sequence.sets <- length( placebo.sequence.sets.list );
  total.num.sequence.sets <- sum( vaccine.num.sequence.sets, placebo.num.sequence.sets );

  dimnames( match.mismatch.table ) <- list( site=1:seq.length, c(paste("match_vax(n1=", vaccine.num.sequence.sets, ")", sep=""), "mismatch_vax",
                                              paste("match_pla(n2=", placebo.num.sequence.sets, ")", sep=""), "mismatch_pla",
                                              paste("match_tot(n=", total.num.sequence.sets, ")", sep=""),"mismatch_tot" ));

  # How many vaccine sequences match the insert at each site? Mismatch at each site?
  match.mismatch.table[ , 1 ] <- sapply( 1:seq.length, function( .seq.i ){ sum( aaPerPosition( vaccine.seq.envir, .seq.i ) == insert.char.vector[ .seq.i ], na.rm = TRUE ) } );
  match.mismatch.table[ , 2 ] <- vaccine.num.sequence.sets - match.mismatch.table[ , 1 ];
  # How many placebo sequences match the insert at each site?
  match.mismatch.table[ , 3 ] <- sapply( 1:seq.length, function( .seq.i ){ sum( aaPerPosition( placebo.seq.envir, .seq.i ) == insert.char.vector[ .seq.i ], na.rm = TRUE ) } );
  match.mismatch.table[ , 4 ] <- placebo.num.sequence.sets - match.mismatch.table[ , 3 ];
  # How many sequences in total match the insert at each site?
  match.mismatch.table[ , 5 ] <- rowSums( match.mismatch.table[ , c( 1, 3 ) ] );
  match.mismatch.table[ , 6 ] <- rowSums( match.mismatch.table[ , c( 2, 4 ) ] );

  return( match.mismatch.table );
} # makeMatchMismatchInsertAAtable (..)

