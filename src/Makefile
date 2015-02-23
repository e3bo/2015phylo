RS = '/usr/bin/Rscript'
MAFFT = '/usr/bin/mafft'
GBLOCKS = '/opt/Gblocks_0.91b/Gblocks'
WD = /root/work
SD = /root/src
DD = /root/data
VPATH = $(SD):$(WD)

aligned.fasta-gb: aligned.fasta
	$(GBLOCKS) $< -t=d -b5=h >gblocks-stdout 2>gblocks-stderr || [ $$? -eq 1 ]

aligned.fasta: pedv-renamed.fasta
	$(MAFFT) --auto $< >$@ 2>mafft-stderr

pedv-renamed.fasta: get-sequences.r pedv-north-america.csv
	$<

pedv-north-america.csv: parse.py pedv-north-america.gb
	$<

pedv-north-america.gb: get-gb.r
	$<
