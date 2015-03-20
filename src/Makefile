RS = /usr/bin/Rscript
MAFFT = /usr/bin/mafft
GBLOCKS = /usr/bin/Gblocks
BEAST = /usr/bin/beast
TREEANN = /usr/bin/treeannotator
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

all: convergence-check.html

convergence-check.html: convergence-check.r
	cp $< . && $(RS) -e "library(knitr); spin(\"$(<F)\")" >$(<F).stdout 2>$(<F).stderr

flagMleReport: mle.r mcc.tree shipment-flows-origins-on-rows-dests-on-columns.csv flagData
	cp $< . && $(RS) -e "library(knitr); spin(\"$(<F)\")" >$(<F).stdout 2>$(<F).stderr && touch $@

shipment-flows-origins-on-rows-dests-on-columns.csv: clean-shipment-matrix.r flagData
	$<

mcc.nh: mcc.tree
	$(RS) -e "library(ape); tree <- read.nexus(\"$<\"); write.tree(tree, file=\"mcc.nh\");"

mcc.tree: pedv.trees
	$(TREEANN) -heights keep -burninTrees 10000 -limit 0.5 $< $@ >treean-stdout 2>treean-stderr

pedv.trees: pedv.xml 
	$(BEAST) -beagle_off $< >beast-stdout 2>beast-stderr

pedv.xml: beast-xml-gen.r flagData sIndel-aligned.fasta-gb nonsIndel-aligned.fasta-gb
	$<

flagData: $(shell find $(DD)/*)
	cp -p $(DD)/* . && touch $@

sIndel-aligned.fasta-gb nonsIndel-aligned.fasta-gb: %-gb: %
	$(GBLOCKS) $< -t=d -b5=h >gblocks-stdout-$* 2>gblocks-stderr-$* || [ $$? -eq 1 ]

sIndel-aligned.fasta nonsIndel-aligned.fasta: flagSequences
	$(MAFFT) --auto $(subst -aligned,,$@) >$@ 2>mafft-stderr-$(subst -aligned.fasta,,$@)

flagSequences: get-sequences.r pedv-north-america.csv
	$< && touch $@

pedv-north-america.csv: parse.py pedv-north-america.gb
	$<

pedv-north-america.gb: get-gb.r
	$<
