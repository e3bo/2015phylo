#!/usr/bin/Rscript

library(ape)
library(rentrez)
library(lubridate)

## sets of accessions based on information in papers
okaCell2014field <- paste0('KM3922', 26:31 + 0.1) # field designation based on Table 1
vlasovaDistinct2014 <- paste0('KJ645', 635:708 + 0.1)
huangOrigin2013 <- paste0('KF46875', 2:4 + 0.1)
marthalerComplete2013 <- 'KF272920.1'
wangNew2014 <- c('KJ399978.1', 'KJ408801.1')
marthalerThird2014 <- 'KM077139.1'
stevensonEmergence2013 <- c('KF452322.1', 'KF452323.1')
chenIsolation2013noPassage <- c('KF650370.1', 'KF650373.1') 
hoangFull2013 <- 'KF804028.1'
wangGenbank <- 'KJ584361.1' # note in Genbank record identifies this as a case

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
foo <- c(vlasovaDistinct2014, okaCell2014field, huangOrigin2013, marthalerComplete2013,
         wangNew2014, marthalerThird2014, stevensonEmergence2013, chenIsolation2013noPassage,
         hoangFull2013, wangGenbank)
feats[foo, 'isFieldSample'] <- TRUE

feats$isPassaged <- FALSE
foo <- c(okaCell2014field, hoangFull2013)
feats[foo, 'isPassaged'] <- TRUE


## from table 1 of Oka et al.
feats['KM392226.1','loc'] <- 'MI'
feats['KM392227.1','loc'] <- 'MI'
feats['KM392228.1', 'loc'] <- 'OH'
feats['KM392229.1', 'loc'] <- 'OH'
feats['KM392230.1', 'loc'] <- 'OH'
feats['KM392231.1', 'loc'] <- 'OH'

## from Genbank record
feats['KF452322.1', 'loc'] <- 'IA'
feats['KF452323.1', 'loc'] <- 'IN'

## from publication
feats['KF650370.1', 'loc'] <- 'IN'
feats['KF650373.1', 'loc'] <- 'IA'
feats[hoangFull2013, 'loc'] <- 'IA'

monyr <- grepl("^\\w{3}-\\d{4}", feats$date)
feats$precision[monyr] <- 30
feats$tipDate[monyr] <- paste0("01-", feats$date[monyr])

dmy <- grepl("^\\d{2}-\\w{3}-\\d{4}$", feats$date)
feats$precision[dmy] <- 1
feats$tipDate[dmy] <- feats$date[dmy]

y <- grepl("^\\d{4}$", feats$date)
foo <- rownames(feats)[y]
last_collection_date_possible <- c("KJ184549.1"=ymd(20131231), "KF468754.1"=ymd(20130726), "KF468753.1"=ymd(20130726))
## These dates based on submission dates in genbank records
end <- last_collection_date_possible[foo]
feats$precision[y] <- as.duration(interval(end=end, start=ymd(20130501)))/ddays(1)
feats$tipDate[y] <- paste0("01-May-", feats$date[y])

query <- paste(rownames(feats), collapse=' ')
res <- entrez_search(db='nucleotide', term=query, retmax=nrow(feats))
seqs <- entrez_fetch(db='nucleotide', rettype='fasta', id=res$ids)
write(seqs, 'pedv.fasta')

dna <- read.dna("pedv.fasta", "fasta")
labs <- names(dna)
labs <- strsplit(labs, split='\\|')
names(dna) <- sapply(labs, '[[', 4)
test <- sapply(dna, length) > 27000

feats$isWholeGenome <- test[rownames(feats)]
feats$tname <- paste(feats$loc, feats$tipDate, feats$precision, rownames(feats),
                     ifelse(feats$isPassaged, 'passaged', 'unpassaged'), sep='_')

test <- feats$isFieldSample %in% TRUE & feats$isWholeGenome
sel <- rownames(feats)[test]

dna <- dna[sel]

names(dna) <- feats[sel, 'tname']

sIndels <- c("IA_20-Oct-2013_1_KJ645649.1_unpassaged",
             "IA_29-Dec-2013_1_KJ645695.1_unpassaged",
             "IA_29-Dec-2013_1_KJ645696.1_unpassaged",
             "IN_08-Jun-2013_1_KJ645635.1_unpassaged",
             "MN_05-Nov-2013_1_KJ645655.1_unpassaged",
             "MN_26-Jun-2013_1_KJ645704.1_unpassaged",
             "OH_29-Jan-2014_1_KJ645702.1_unpassaged",
             "OH_15-Jan-2014_1_KJ399978.1_unpassaged")

write.dna(dna[sIndels], file='sIndels-renamed', format='fasta')
test <- !names(dna) %in% sIndels
write.dna(dna[test], file='nonsIndels.fasta', format='fasta')
