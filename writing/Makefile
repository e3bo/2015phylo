RS = '/usr/bin/Rscript'

main-text.pdf: main-text.tex
	latexmk -pdf $<

main-text.tex: main-text.Rnw
	$(RS) -e "knitr::knit(\"main-text.Rnw\")"
