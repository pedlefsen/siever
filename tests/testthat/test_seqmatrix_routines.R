library("testthat")
library("siever")
library("seqinr")

context("Testing seq matrix reading routines")
test_that(
   "Alignment fasta files with differing length alignments generate errors",
   {expect_that(siever:::getSeqMtx("control_A-bad.fasta"),throws_error())}
)
test_that(
   "Alignment fasta files with differing length alignments can generate warnings and recover.",
   {expect_that(siever:::getSeqMtx("control_A-bad.fasta",warn.and.recover=T),gives_warning())}
)

test_that(
   "Matrix filtering works as expected",
   {
      mtx <- siever:::getSeqMtx("control_A.fasta")
      expect_that(length(grep("\\*",mtx))==0,is_false())
      filtmtx <- siever:::filterAcceptable(mtx)
      expect_that(length(grep("\\*",filtmtx))==0,is_true())   
   }
)
      
