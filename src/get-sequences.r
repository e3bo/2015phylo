#!/usr/bin/Rscript

library(ape)
library('rentrez')
library(lubridate)

## sets of acceccsions based on information in papers
okaCell2014field <- paste0('KM3922', 26:31 + 0.1) # field designation based on Table 1
vlasovaDistinct2014 <- paste0('KJ645', 635:708 + 0.1)
huangOrigin2013 <- paste0('KF46875', 2:4 + 0.1)
marthalerComplete2013 <- 'KF272920.1'
wangNew2014 <- 'KJ399978.1'
## identification of 3-rd strain

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

feats$isFieldSample <- NA
feats[vlasovaDistinct2014, 'isFieldSample'] <- TRUE
feats[okaCell2014field, 'isFieldSample'] <- TRUE
feats[huangOrigin2013, 'isFieldSample'] <- TRUE
feats[marthalerComplete2013, 'isFieldSample'] <- TRUE
feats[wangNew2014, 'isFieldSample'] <- TRUE
feats['KM392226.1','loc'] <- 'MI'
feats['KM392227.1','loc'] <- 'MI'
feats['KM392228.1', 'loc'] <- 'OH'
feats['KM392229.1', 'loc'] <- 'OH'
feats['KM392230.1', 'loc'] <- 'OH'
feats['KM392231.1', 'loc'] <- 'OH'

monyr <- grepl("^\\w{3}-\\d{4}", feats$date)
feats$precision[monyr] <- 30
feats$tipDate[monyr] <- paste0("01-", feats$date[monyr])

dmy <- grepl("^\\d{2}-\\w{3}-\\d{4}$", feats$date)
feats$precision[dmy] <- 1
feats$tipDate[dmy] <- feats$date[dmy]

y <- grepl("^\\d{4}$", feats$date)
foo <- feats$accession[y]
last_collection_date_possible <- c("KJ184549.1"=ymd(20131231), "KF468754.1"=ymd(20130726), "KF468753.1"=ymd(20130726))
## These dates based on submission dates in genbank records
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
