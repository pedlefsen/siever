# THIS FUNCTION TESTS THE gwj_analyze WITH AND WITHOUT WRAPPER AGAINST THE CURRENT IMPLEMENTATION OF T-TEST
# CURRENTLY SETUP = 0 (DEFAULT) AND SETUP = 1 ARE ONLY FUNCTIONAL
# SET THE TWO PATH ACCORDING TO WHERE THE CORRESPONDING FILES ARE LOCALLY ON YOUR MACHINE
library("testthat")
library("siever")
source("gwj_methods.R")

test_gwjAnalyze <- function (setup=0){

  vac.file <- system.file("extdata","control_A.fasta",package="siever")
  plc.file <- system.file("extdata","control_B.fasta",package="siever")
  ins.file <- system.file("extdata","control_HXB2.fasta",package="siever")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEFAULT SETUP FOR THE PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vaccine.seq.envir <- readFastaFile(.file.name = paste(vac.file), be.verbose = FALSE)
  placebo.seq.envir <- readFastaFile(.file.name = paste(plc.file), be.verbose = FALSE)
  insert.seq.envir <- readFastaFile(.file.name = paste(ins.file), be.verbose = FALSE)
  weight.matrix = NULL 
  site.sets.list = NULL
  use.f.test = FALSE 
  return.t.test.result = FALSE
  weights.across.sites.in.a.set.init = 0 
  weights.across.sites.in.a.set.accumulation.fn = sum
  vaccine.sequence.sets.list = NULL
  placebo.sequence.sets.list = NULL 
  mimic.smmb = FALSE 
  weights.across.sequences.in.a.set.accumulation.fn = if( mimic.smmb ) { sum } else { mean }
  if( is.null( weight.matrix ) ) { 
    acceptable.chars = c( "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-" )
  } else{ acceptable.chars = setdiff( rownames( weight.matrix ), 'X' ) }
  instead.return.weights = FALSE
  
  insert.char.vector <- strsplit(as.list(insert.seq.envir)[[1]],"")[[1]]
  include.sites <- NULL
  choose.sites.func <- choose.sites.func.arglist <- NULL
  KL <- SKL <- FALSE
  
  # NOTE: WARNING ALL SEQS MUST HAVE SAME LENGTH IN FASTA FILES, BUT HERE THE INPUT DOESN'T FOR A COUPLE OF CASES AND THOSE
  #       DISCARDED. A WARNING IS GENERATED TO NOTIFY OF THE DISCARDING. THIS IS SOLELY DUE TO BAD INPUT AND HAS NO IMPACT 
  #       ON TESTING. ALL FUNCTIONS ARE TESTED WITH THE "CORRECTED" FASTA FILES AND THEREFORE WARNING IS SUPRESSED
  suppressWarnings(vaccine.seq.mtx <- siever:::convert.env.to.matrix(seq.env=vaccine.seq.envir)) 
  suppressWarnings(placebo.seq.mtx <- siever:::convert.env.to.matrix(seq.env=placebo.seq.envir)) 
  
  # Note: using newer, seqinr-dependent routines to create data matrix
  vsm <- siever:::filterAcceptable(siever:::getSeqMtx(vac.file,prefix=""),replace.with="-")
  psm <- siever:::filterAcceptable(siever:::getSeqMtx(plc.file,prefix=""),replace.with="-")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SETUP 1 TESTS THE FUNCTIONALITY OF SUMMING OVER SITES AND OBSERVATIONS SIMULTANEOUSLY
  # NOTE THE LISTS ARE ARBITARILY MADE JUST FOR THE SAKE OF TESTING
  if(setup == 1){
    site.sets.list <- list(c(10:11), c(12:30), c(100:110)) 
    vac.obs.id <- ls(envir = vaccine.seq.envir)
    vaccine.sequence.sets.list <- list(vac.obs.id[50:60], vac.obs.id[200:400], vac.obs.id[500:510])
    plc.obs.id <- ls(envir = placebo.seq.envir)
    placebo.sequence.sets.list <- list(plc.obs.id[50:60], plc.obs.id[200:400], plc.obs.id[420:430])
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mimic OPTION WITHOUT ACCUMMULATION
  if(setup == 2){
    mimic.smmb <- T
    weights.across.sequences.in.a.set.accumulation.fn <- sum 
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mimic OPTION WITH ACCUMMULATION 
  if(setup == 3){
    site.sets.list <- list(c(10:11), c(12:30), c(100:110)) 
    vac.obs.id <- ls(envir = vaccine.seq.envir)
    vaccine.sequence.sets.list <- list(vac.obs.id[50:60], vac.obs.id[200:400], vac.obs.id[500:510])
    plc.obs.id <- ls(envir = placebo.seq.envir)
    placebo.sequence.sets.list <- list(plc.obs.id[50:60], plc.obs.id[200:400], plc.obs.id[420:430])
    mimic.smmb <- T
    weights.across.sequences.in.a.set.accumulation.fn <- sum
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP 4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(setup == 4){
    KL <- TRUE
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SETUP 5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(setup == 5){
    SKL <- TRUE
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  vaccine.indices <- 1:nrow(vaccine.seq.mtx)
  dmtx <- rbind(vsm,psm)
  
  stat.type <- "tStat"
  if(mimic.smmb) stat.type <- "mimic.smmb" 
  if(KL) stat.type <- "KL"
  if(SKL) stat.type <- "SKL"
  
  
  wOpt <-
    weightopt(weight.mtx = weight.matrix, stat.type= stat.type, insert.char.vector = insert.char.vector,
              site.sets.list= site.sets.list, obs.vSeq.sets.list= vaccine.sequence.sets.list,
              obs.pSeq.sets.list= placebo.sequence.sets.list, acceptable.chars = acceptable.chars,
              within.sites.acc.fn = weights.across.sites.in.a.set.accumulation.fn,
              within.obs.acc.fn = weights.across.sequences.in.a.set.accumulation.fn,
              choose.sites.func = choose.sites.func, choose.sites.func.arglist = choose.sites.func.arglist)
  
  if(wOpt$stat.type == "KL"){
    ptm <- proc.time()
    new.ks <- GWJsieve(data.mtx=dmtx ,vaccine.indices=vaccine.indices,wOpt = wOpt,sw.test=T)
    analyze.time <- proc.time()- ptm
    
    ptm <- proc.time()
    current.ks <- klStat(vaccine.seq.envir= vaccine.seq.envir , placebo.seq.envir = placebo.seq.envir )
    current.time <- proc.time() - ptm
    
    cat("Set up was ", setup, "\n")
    cat("Current KL stat...............time = ", round(current.time[3],2), "s\n", "head: " ,head(current.ks),"\n")
    cat("new KL stat...................time = ", round(analyze.time[3],2), "s\n", "head: ", head(new.ks),  ifelse(isTRUE(all.equal(current.ks, new.ks,check.attributes=F)),"All matched\n", "Didn't match\n"))
    
    expect_that(new.ks,is_equivalent_to(current.ks))
    return()
  }
  
  if(wOpt$stat.type == "SKL"){
    ptm <- proc.time()
    new.sks <- GWJsieve(data.mtx=dmtx ,vaccine.indices=vaccine.indices,wOpt = wOpt,sw.test=T)
    analyze.time <- proc.time()- ptm
    
    ptm <- proc.time()
    current.sks <- sklStat(vaccine.seq.envir= vaccine.seq.envir , placebo.seq.envir = placebo.seq.envir )
    current.time <- proc.time() - ptm
    
    cat("Set up was ", setup, "\n")
    cat("Current KLS stat...............time = ", round(current.time[3],2), "s\n", "head: " ,head(current.sks),"\n")
    cat("new KLS stat...................time = ", round(analyze.time[3],2), "s\n", "head: ", head(new.sks),  ifelse(isTRUE(all.equal(current.sks, new.sks,check.attributes=F)),"All matched\n", "Didn't match\n"))
    
    expect_that(new.sks,is_equivalent_to(current.sks))
    return()
  }
  
  
  ptm <- proc.time()
  
  tst.alt <- 
    GWJsieve(data.mtx=dmtx, vaccine.indices=vaccine.indices,wOpt=wOpt,sw.test=T)
  
  analyze.time <- proc.time() - ptm
  
  #ptm <- proc.time()
  #tst.alt.5 <- 
  #  GWJsieve(data.mtx=dmtx, vaccine.indices=vaccine.indices, perm.num=5, wOpt=wOpt)
  #permute.time <- proc.time() - ptm
  
  ptm <- proc.time()
  wrapper.analyze <-
    GWJsieve.wrapper(
      vaccine.seq.envir= vaccine.seq.envir, placebo.seq.envir = placebo.seq.envir,
      insert.seq.envir = insert.seq.envir, weight.matrix = weight.matrix,
      site.sets.list = site.sets.list, use.f.test = use.f.test,
      return.t.test.result = return.t.test.result,
      weights.across.sites.in.a.set.init = weights.across.sites.in.a.set.init,
      weights.across.sites.in.a.set.accumulation.fn = weights.across.sites.in.a.set.accumulation.fn,
      vaccine.sequence.sets.list = vaccine.sequence.sets.list,
      placebo.sequence.sets.list = placebo.sequence.sets.list,
      mimic.smmb = mimic.smmb,
      weights.across.sequences.in.a.set.accumulation.fn = weights.across.sequences.in.a.set.accumulation.fn,
      acceptable.chars = acceptable.chars, instead.return.weights = instead.return.weights
    )
  wrapper.time <- proc.time() - ptm
  
  ptm <- proc.time()
  current.tStat <-
    tStat(
      vaccine.seq.envir= vaccine.seq.envir, placebo.seq.envir = placebo.seq.envir,
      insert.seq.envir = insert.seq.envir, weight.matrix = weight.matrix,
      site.sets.list = site.sets.list, use.f.test = use.f.test,
      return.t.test.result = return.t.test.result,
      weights.across.sites.in.a.set.init = weights.across.sites.in.a.set.init,
      weights.across.sites.in.a.set.accumulation.fn = weights.across.sites.in.a.set.accumulation.fn,
      vaccine.sequence.sets.list = vaccine.sequence.sets.list,
      placebo.sequence.sets.list = placebo.sequence.sets.list,
      mimic.smmb = mimic.smmb,
      weights.across.sequences.in.a.set.accumulation.fn = weights.across.sequences.in.a.set.accumulation.fn,
      acceptable.chars = acceptable.chars, instead.return.weights = instead.return.weights
    )
  current.time <- proc.time() - ptm
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRINTING THE TEST RESULTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cat("Set up was ", setup, "\n")
  cat("Current tStat...............time = ", round(current.time[3],2), "s\n", "head: " ,head(current.tStat),"\n")
  cat("Wrapper result..............time = ", round(wrapper.time[3],2), "s\n", "head: ", head(wrapper.analyze),ifelse(isTRUE(all.equal(current.tStat, wrapper.analyze,check.attributes=F)),"All matched\n", "Didn't match\n"))
  cat("new tStat...................time = ", round(analyze.time[3],2), "s\n", "head: ", head(tst.alt),  ifelse(isTRUE(all.equal(current.tStat, tst.alt,check.attributes=F)),"All matched\n", "Didn't match\n"))
  #cat("new permute tStat:..........time = ", round(permute.time[3],2), "s\n", "1st 6 (or less) elements: ","\n")
  
  #print(tst.alt.5[1:5,1:min(6,ncol(tst.alt.5))])
  
  expect_that(tst.alt,is_equivalent_to(current.tStat))
  expect_that(wrapper.analyze,is_equivalent_to(current.tStat))
  
#  expect_that(current.time[3]>wrapper.time[3],is_true())
#  expect_that(current.time[3]>analyze.time[3],is_true())
    
}

#RUNNING THE TEST
context("Does the new tStat method return the same result as the old one?")
#default
test_that(
   "Original tStat yields the same result as the new one, with setup 0",
   {test_gwjAnalyze()}
) 
# WITH ACCUMULATION ACROSS SITES AND OBSERVATIONS
test_that(
   "Original tStat yields the same result as the new one with setup 1",
   {test_gwjAnalyze(setup = 1)}
)

# SMMB WEIGHT OPTION
test_that(
  "Original tStat yields the same result as the new one with setup 2",
{test_gwjAnalyze(setup = 2)}
)

# SMMB WEIGHT OPTION WITH ACCUMULATION
test_that(
  "Original tStat yields the same result as the new one with setup 3",
{test_gwjAnalyze(setup = 3)}
)

# KL Option
test_that(
  "Original tStat yields the same result as the new one with setup 4",
{test_gwjAnalyze(setup = 4)}
)

# SKL Option
test_that(
  "Original tStat yields the same result as the new one with setup 5",
{test_gwjAnalyze(setup = 5)}
)