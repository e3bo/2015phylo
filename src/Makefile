RS = '/usr/bin/Rscript'
MAFFT = '/usr/bin/mafft'
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

all: aligned.fasta

aligned.fasta: pedv-renamed.fasta
	$(MAFFT) --auto $< >$@

pedv-renamed.fasta: get-sequences.r pedv-north-america.csv
	$<

pedv-north-america.csv: parse.py pedv-north-america.gb
	$<

pedv-north-america.gb: get-gb.r
	$<
