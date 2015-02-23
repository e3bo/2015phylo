#!/usr/bin/Rscript

library(ape)
library('rentrez')
library(lubridate)

feats <- read.csv('pedv-north-america.csv', stringsAsFactors=FALSE)
foo <- strsplit(feats$country, split=':[ ]*')
tmpf <- function(x) x[[length(x)]]
foo <- sapply(foo, tmpf)
foo[foo=='NorthCarolina'] <- 'North Carolina'
foo[foo=='Tennesse'] <- 'Tennessee'
foo[foo=='CO'] <- 'Colorado'
stopifnot(foo %in% c(state.name, 'Canada', 'USA', 'Mexico'))
ind <- match(foo, state.name)
feats$abb <- state.abb[ind]
feats$loc <- ifelse(is.na(feats$abb), feats$country, feats$abb)

monyr <- grepl("^\\w{3}-\\d{4}", feats$date)
feats$precision[monyr] <- 30
feats$tipDate[monyr] <- paste0("01-", feats$date[monyr])

dmy <- grepl("^\\d{2}-\\w{3}-\\d{4}$", feats$date)
feats$precision[dmy] <- 1
feats$tipDate[dmy] <- feats$date[dmy]

y <- grepl("^\\d{4}$", feats$date)
foo <- feats$accession[y]
last_collection_date_possible <- c("KJ184549.1"=ymd(20131231), "KF468754.1"=ymd(20130902), "KF468753.1"=ymd(20130902))
## These dates based on publication history of any papers linked in genbank records
end <- last_collection_date_possible[foo]
feats$precision[y] <- as.duration(interval(end=end, start=ymd(20130501)))/ddays(1)
feats$tipDate[y] <- paste0("01-May-", feats$date[y])

query <- paste(feats$accession, collapse=' ')
res <- entrez_search(db='nucleotide', term=query, retmax=100)
seqs <- entrez_fetch(db='nucleotide', rettype='fasta', id=res$ids)
write(seqs, 'pedv.fasta')

dna <- read.dna("pedv.fasta", "fasta")
labs <- names(dna)
labs <- strsplit(labs, split='\\|')
acsn <- sapply(labs, '[[', 4)
feats$tname <- paste(feats$loc, feats$tipDate, feats$precision, feats$accession, sep='_')
ind <- match(acsn, feats$accession)
test <- sapply(dna, length) > 27000
names(dna) <- feats$tname[ind]
write.dna(dna[test], file='pedv-renamed.fasta', format='fasta')
