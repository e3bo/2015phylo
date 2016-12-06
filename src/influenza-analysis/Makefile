all: influenza-c5.RData

influenza-c5.RData: fit-model-selected-from-simulated-data.R influenza-c4.RData
	$< >fit-mod-sim.stdout 2>fit-mod-sim.stderr

influenza-c4.RData: stability-selection-on-simulated-data.R influenza-c3.RData
	$< >stab-sel-sim.stdout 2>stab-sel-sim.stderr

influenza-c3.RData: fit-model-selected-from-real-data.R influenza-c2.RData
	$< >fit-mod.stdout 2>fit-mod.stderr

influenza-c2.RData: stability-selection-on-real-data.R influenza-c1.RData
	$< >stab-sel.stdout 2>stab-sel.stderr

influenza-c1.RData: estimate-tree-from-alignment.R
	$< >est-trees.stdout 2>est-trees.stderr

clean:
	rm influenza-c[12345].RData Rplots.pdf *.stdout *.stderr
