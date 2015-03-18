#!/usr/bin/Rscript

library(ape)
library(lubridate)
library(XML)

xmlTemplatePath <- file.path("beast-template.xml")
outputBasename <- "pedv"
alignmentPath <- file.path("aligned.fasta-gb")
chainLength <- 1e8
logEvery <- 1e4

aln <- read.dna(alignmentPath, format='fasta')
taxaNames <- rownames(aln)

xmlTree <- xmlTreeParse(xmlTemplatePath, getDTD=F)
root <- xmlRoot(xmlTree)

days2years <- function(x) as.integer(x) * ddays(1) / dyears(1)

tmpf <- function(x) {
    fields <- strsplit(x, "_")[[1]]
    tipdate <- decimal_date(dmy(fields[2]))
    windowLength <- days2years(fields[3])
    xmlNode("taxon", attrs=c(id=x),
            xmlNode("date", attrs=c(value=as.character(tipdate),
                                direction="forwards", units="years",
                                precision=as.character(windowLength)))
            )
}
taxaNodeList <- lapply(taxaNames, tmpf)
taxaNode <- xmlNode("taxa", attrs=c(id="taxa"))
root[["taxa"]] <- addChildren(taxaNode, kids=taxaNodeList)

tmpf <- function(x) {
    seq <- toupper(paste(as.character(x)[1,], collapse=''))
    idref <- labels(x)
    xmlNode("sequence", xmlNode("taxon", attrs=c(idref=idref)), xmlTextNode(seq))
}

sequenceNodeList <- lapply(aln, tmpf)
alignmentNode <- xmlNode("alignment", attrs=c(id="alignment", dataType="nucleotide"))
root[["alignment"]] <- addChildren(alignmentNode, kids=sequenceNodeList)

ind <- which(names(xmlChildren(root[["mcmc"]])) == "log")

for(i in ind){
    x <- root[["mcmc"]][[i]]
    attrs <- xmlAttrs(x)
    if("fileName" %in% names(attrs)){
        xmlAttrs(root[["mcmc"]][[i]])['fileName'] <- paste0(outputBasename, '.log')
        xmlAttrs(root[["mcmc"]][[i]])['logEvery'] <- as.character(as.integer(logEvery))
    }
}

xmlAttrs(root[["mcmc"]][["logTree"]])["fileName"] <- paste0(outputBasename, '.trees')
xmlAttrs(root[["mcmc"]])["operatorAnalysis"] <- paste0(outputBasename, '.ops')
xmlAttrs(root[["mcmc"]])["chainLength"] <- as.character(as.integer(chainLength))

sink(paste0(outputBasename, '.xml'))
root
sink()      

q('no')

## Rest is all about allowing variable tip dates, which seems to lead
## to negativel coalescent for some reason

#'  | reference | doi | gene | estimate | HPD |
#'  |-----------------------------------------|
#'  | Vijgen et al. 2006 | 10.1128/JVI.02675-05 | PHEV, BCoV, HCoV-OC43 | spike | 6.1E-4 | 2.1E-4 | 10E-4 |
#'  |                    |                      |                       | nucleocapsid | 3.6e-4 | 1.1e-4 | 6.3e-4 |
#'  | Vijgen et al. 2005 | 10.1128/JVI.79.3.1595-1604.2005 | BCoV, HCoV-OC43 | spike | 4.3e-4 | 2.7e-4 | 6e-4 |
#'  | Sanchez et al. 1992 |                     |  TGEV |  not in abstract, probably spike |  7e-4 | 5e-4 | 9e-4 |


tmpf <- function(x){
    parmid <- paste0(x, ".height")
    xmlNode("leafHeight", attrs=c(taxon=x), xmlNode("parameter", attrs=c(id=parmid)))
}

leafHeightList <- lapply(taxaNames, tmpf)
root[['treeModel']] <- addChildren(root[['treeModel']], kids=leafHeightList)

tmpf <- function(x){
    parmid <- paste0(x, ".height")
    fields <- strsplit(x, "_")[[1]]
    windowSize <- days2years(fields[3])*0.1
    xmlNode("randomWalkOperator", attrs=c(windowSize=as.character(windowSize),
                                      weight="1"),
            xmlNode("parameter", attrs=c(idref=parmid)))
}

rwoList <- lapply(taxaNames, tmpf)
root[["operators"]] <- addChildren(root[["operators"]], kids=rwoList)

tmpf <- function(x){
    parmid <- paste0(x, ".height")
    xmlNode("column", attrs=c(label=parmid, sf=3, width=6),
      xmlNode("parameter", attrs=c(idref=parmid)))
}

lhpList <- lapply(taxaNames[c(83:86)], tmpf)
root[["mcmc"]][["log"]] <- addChildren(root[["mcmc"]][["log"]], kids=lhpList)

