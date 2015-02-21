#!/usr/bin/Rscript

library('rentrez')

search <- entrez_search(db="nucleotide", term="Porcine epidemic diarrhea virus[Organism] AND (country=USA OR country=Mexico OR country=Canada)", retmax=150)
nuc <- entrez_fetch(db="nucleotide", id=search$ids, rettype="gb")
write(nuc, file='pedv-north-america.gb')
