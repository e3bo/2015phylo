RS = '/usr/bin/Rscript'

main-text.odt: main-text.md pedv.bib
	pandoc --filter pandoc-citeproc -t odt -o main-text.odt $<

main-text.md: main-text.Rmd
	$(RS) -e "knitr::knit(\"main-text.Rmd\")"
