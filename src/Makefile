SHELL = bash
RS = /usr/bin/Rscript
MAFFT = /usr/bin/mafft
GBLOCKS = /usr/bin/Gblocks
BEAST = /usr/bin/beast
TREEANN = /usr/bin/treeannotator
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

all: flagMleReport posteriors-and-priors.pdf

posteriors-and-priors.pdf: plot-posteriors-and-priors.r beast/priors/pedv-priors.log beast/run1/pedv.log beast/run2/pedv.log beast/run3/pedv.log beast/run4/pedv.log nonsIndel-aligned.fasta-gb sIndel-aligned.fasta-gb
	$<

flagMleReport: mle.r convergence-check.html shipment-flows-origins-on-rows-dests-on-columns.csv flagData state-hogBalanceSheetDec2011Dec2012.csv
	cp $< . && $(RS) -e "library(knitr); spin(\"$(<F)\")" >$(<F).stdout 2>$(<F).stderr && touch $@

convergence-check.html: convergence-check.r beast/run1/pedv.log beast/run2/pedv.log beast/run3/pedv.log beast/run4/pedv.log
	cp $< . && $(RS) -e "library(knitr); spin(\"$(<F)\")" >$(<F).stdout 2>$(<F).stderr

shipment-flows-origins-on-rows-dests-on-columns.csv: clean-shipment-matrix.r flagData
	$<

mcc.nh: mcc.tree
	$(RS) -e "library(ape); tree <- read.nexus(\"$<\"); write.tree(tree, file=\"mcc.nh\");"

mcc.tree: run1/pedv.trees
	$(TREEANN) -heights keep -burninTrees 10000 -limit 0.5 $< $@ >treean-stdout 2>treean-stderr

beast/priors/pedv-priors.log: flagData
	bd=$(subst pedv-priors.log,,$@) && install -d $$bd && cp pedv-priors.xml $$bd \
	&& $(BEAST) -working -beagle_off $${bd}pedv-priors.xml >$${bd}beast-stdout 2>$${bd}beast-stderr &\
	while [ ! -f $(subst pedv-priors.log,,$@)pedv-priors.ops ]; do sleep 60; done; kill $$(jobs -p)

beast/run1/pedv.log beast/run2/pedv.log beast/run3/pedv.log beast/run4/pedv.log: sample-initial-values-from-priors.r pedv.xml
	bd=$(subst pedv.log,,$@) && install -d $$bd && $< > $${bd}pedv-reinitialized.xml \
	&& $(BEAST) -working -beagle_off $${bd}pedv-reinitialized.xml >$${bd}beast-stdout 2>$${bd}beast-stderr &\
	while [ ! -f $(subst pedv.log,,$@)pedv.ops ]; do sleep 60; done; kill $$(jobs -p)

pedv.xml: flagData sIndel-aligned.fasta-gb nonsIndel-aligned.fasta-gb
	md5sum -c alignment-md5sums && touch $@

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

state-hogBalanceSheetDec2011Dec2012.csv: get-balance-sheet-2012-csv.sh
	$<

clean:
	rm convergence-check.r mle.r
