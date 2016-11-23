#!/usr/bin/Rscript

library(ape)
library(lubridate)

library(XML)

xmlTemplatePath <- file.path("../data/origins.xml")

xmlTree <- XML::xmlParse(xmlTemplatePath)
treel <- XML::xmlToList(xmlTree)

align_els <- c(H1N1=11, H1N2=14)

get_dnabin <- function(align_id){
    aln <- treel[[align_id]]
    aln <- aln[-length(aln)]
    seqs <- lapply(aln, "[[", "text")
    names(seqs) <- sapply(aln, "[[", "taxon")
    stripws <- function(x){
      gsub("\t|\n", "", x)
    }
    seqs <- lapply(seqs, stripws)

    ncol <- as.integer(sapply(seqs, nchar))
    stopifnot(all(ncol == ncol[1]))
    ncol <- ncol[1]
    nseq <- length(seqs)
    tf <- tempfile()
    cat(nseq, " " , ncol, "\n", sep="", file=tf)
    for(i in seq(1, nseq)){
      cat(names(seqs)[i], "\t", seqs[[i]], "\n", file = tf, append=TRUE)
    }
    read.dna(tf, format = "sequential")
}

almnts <- lapply(align_els, get_dnabin)

saveRDS(almnts, file="../data/swine-influenza-alignments-dnabin.rds")
q("no")
