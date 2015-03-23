#!/usr/bin/Rscript

library(XML)

tree <- xmlParse("pedv.xml", getDTD=F)

rval <- rexp(n=1,rate=2)
xpathApply(tree, "/beast/siteModel/gammaShape/parameter", 'xmlAttrs<-', value=c(value=rval))

rval <- exp(rnorm(n=1, mean=1, sd=1.25))
xpathApply(tree, "/beast/HKYModel/kappa/parameter", 'xmlAttrs<-', value=c(value=rval))

rval <- runif(4)
rval <- rval/sum(rval)
rval <- paste(rval, collapse=' ')
xpathApply(tree, "/beast/HKYModel/frequencies/frequencyModel/frequencies/parameter", 'xmlAttrs<-', value=c(value=rval))

rval <- rnorm(n=1, mean=7e-4, sd=1e-4)
while(rval <= 0){
    rval <- rnorm(n=1, mean=7e-4, sd=1e-4)
}    
xpathApply(tree, "/beast/strictClockBranchRates/rate/parameter", 'xmlAttrs<-', value=c(value=rval))

saveXML(xmlRoot(tree), file='pedv-reinitialized.xml')
